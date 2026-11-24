#!/usr/bin/env python3
"""cloudprofile_pv.py -- ParaView Python View: cloud / ice densities against height, labelled T [K] (left) and p [bar] (right).

One panel per latitude, read from a zonal slice <Planet>_zonal_<k>_<it>.vtk (the model is zonally symmetric, so a column
of the slice is the column at that latitude):
  x axis   condensate density [g/m3], log: liquid clouds SOLID, ice DASHED, one colour per species
           (H2O, NH3, CH4, H2S where the slice has it), NH4SH (the chemistry's solid product) dotted
  y axis   the model's level height (uniform); the LEFT ticks print the temperature [K] at those levels, the RIGHT ticks
           the static pressure [bar] at the same levels -- labels, not separate scales, so a T that does not fall
           monotonically (a stratosphere) stays readable
  shading  the model's mixed-phase band of each species, t_00 .. t_0 (all four giants: H2O 210.15..273.15 K, NH3
           140..195.5 K, CH4 67.36..90.69 K; Uranus / Neptune also H2S 140..187.65 K; the <TAG>_T00_LEGACY values
           are not drawn): below t_00 the condensate is all ice, above t_0 all liquid

Atreya overlay (2026-11-24, --atreya): Atreya's equilibrium cloud models (ECCM), digitised by
literature/atreya_digitize/digitize.py into literature/atreya_eccm.csv, drawn behind the model curves in the species' colour:
  shaded + base tick  placed BY PRESSURE: the deck's p from the figure's right axis -> the model level with that p
  dashed line         the same curve placed BY TEMPERATURE: Atreya's T -> the level where the model's troposphere has that T
Densities x1000 (Atreya plots g/L = kg/m3, the model g/m3); they are his upper limits. Cases: Jupiter 1x (A&W 2005 Fig. 2),
Saturn 5x (Fig. 3), Uranus solar (Atreya 2020 Fig. 5), Neptune 30x (A&W Fig. 5), Uranus / Neptune SN6_N0.001 / SN6_N0.01 /
SN6_N0.1 / SN16_N5 (Atreya 2020 Fig. 6). Default per planet: Jupiter 1x, Saturn 5x, Uranus and Neptune SN16_N5 (the set-up
S/N is ~10x solar, Fig. 6's regime); --atreya none switches it off. The csv is looked up at ../../literature/atreya_eccm.csv
from this file, or in ATREYA_ECCM.

Works for all four giants (2026-11-21): the planet comes from the file name; Saturn writes T in K/10, Jupiter /
Neptune / Uranus in degC; NH4SH is written x1e9 (ug/m3) everywhere. The x range follows the largest density
(six decades). Default latitudes 0, 45, 75 on every planet; --lat for anything else. The SAME file is in every
giant's tools/.

Command line (system python3, which has ParaView 5.10 + matplotlib; pvpython has no matplotlib here):
  python3 -I tools/cloudprofile_pv.py <zonal.vtk> [--lat 0,45,75] [--out file.png] [--size 1800x900] [--atreya <case>|none]
  (output default: <zonal>_cloudprofile.png next to the slice)
Inside the ParaView GUI: View > Python Shell > Run Script, choose this file; it asks nothing and plots the active
source if that is a zonal slice, else opens the file in CLOUDPROFILE_FILE (environment) -- or edit DEFAULT_FILE below.
The view is a normal ParaView Python View: File > Save Screenshot works, and the script can be edited in the view's
properties (Script).
"""
import os, sys

DEFAULT_FILE = ""          # used inside the GUI when no zonal slice is the active source

# ----------------------------------------------------------------------------------------------- the Python View script
VIEW_SCRIPT = r'''
LATS = @LATS@
PLANET = @PLANET@
TITLE = @TITLE@
ATREYA = @ATREYA@          # [(species, phase, [T], [p], [rho g/m3]), ...] or [] ; ATREYA_LABEL names the case / figure
ATREYA_LABEL = @ATREYA_LABEL@
FS_TICK, FS_LABEL, FS_TITLE, FS_LEGEND, FS_BAND, FS_SUPTITLE = 13, 15, 16, 14, 12, 16   # text sizes (2026-11-21: larger, user)

def setup_data(view):
    from paraview.vtk import vtkDataObject
    for n in range(view.GetNumberOfVisibleDataObjects()):
        ds = view.GetVisibleDataObjectForSetup(n)
        pd = ds.GetPointData()
        for a in range(pd.GetNumberOfArrays()):
            view.SetAttributeArrayStatus(n, vtkDataObject.POINT, pd.GetArrayName(a), 1)

def render(view, width, height):
    import numpy as np
    from paraview import python_view
    from vtk.util.numpy_support import vtk_to_numpy
    fig = python_view.matplotlib_figure(width, height)
    if view.GetNumberOfVisibleDataObjects() < 1:
        return python_view.figure_to_image(fig)
    ds = view.GetVisibleDataObjectForRendering(0)
    pd = ds.GetPointData()
    npt = ds.GetNumberOfPoints()
    dims = ds.GetDimensions() if hasattr(ds, "GetDimensions") else (0, 0, 0)
    nj, ni = dims[0], dims[1]                            # zonal slice: DIMENSIONS jm im 1, point index i*jm + j
    if nj * ni != npt:                                   # the view's delivered copy can come without its extent
        hh = pd.GetArray("height")
        ni = len(np.unique(vtk_to_numpy(hh))) if hh is not None else 41
        nj = npt // ni
    def arr(n):
        a = pd.GetArray(n)
        return None if a is None else vtk_to_numpy(a).astype(float).reshape(ni, nj)
    T = arr("Temperature"); P = arr("PressureStat"); H = arr("height")
    T = T * 10.0 if PLANET == "Saturn" else T + 273.15  # Saturn writes K/10, the other giants degC
    if H is None: H = np.tile(np.arange(ni, dtype=float)[:, None], (1, nj))
    sp = [  # (label, liquid array, ice array, colour, factor to g/m3)
        ("H2O",   "H2OCloud", "H2OIce", "#1f77b4", 1.0),
        ("NH3",   "NH3Cloud", "NH3Ice", "#2ca02c", 1.0),
        ("CH4",   "CH4Cloud", "CH4Ice", "#d62728", 1.0),
        ("H2S",   "H2SCloud", "H2SIce", "#9467bd", 1.0)]
    base_b = [("H2O", 210.15, 273.15, "#1f77b4"), ("NH3", 140.0, 195.5, "#2ca02c"), ("CH4", 67.36, 90.69, "#d62728")]
    h2s_b = [("H2S", 140.0, 187.65, "#9467bd")]
    bands = {"Saturn": base_b, "Jupiter": base_b, "Neptune": base_b + h2s_b, "Uranus": base_b + h2s_b}.get(PLANET, base_b)
    dmax = 0.0                                           # the x range follows the largest condensate density
    for (nm, liq, ice, c, f) in sp:
        for n_ in (liq, ice):
            a_ = arr(n_)
            if a_ is not None and np.isfinite(a_).any(): dmax = max(dmax, float(np.nanmax(a_)) * f)
    S_ = arr("NH4SH")
    if S_ is not None and np.isfinite(S_).any(): dmax = max(dmax, float(np.nanmax(S_)) * 1e-6)
    for d_ in ATREYA: dmax = max(dmax, max(d_[4]))
    acol = {"H2O": "#1f77b4", "NH3": "#2ca02c", "CH4": "#d62728", "H2S": "#9467bd", "NH4SH": "#8c564b"}
    lats = [l for l in LATS if -90 <= l <= 90]
    axs = fig.subplots(1, len(lats), sharey=False, squeeze=False)[0]   # own T labels per latitude
    xmax = 10.0 ** np.ceil(np.log10(dmax)) if dmax > 0 else 1e2; xmin = xmax * 1e-6
    for ax, lat in zip(axs, lats):
        j = int(round((90.0 - lat) * (nj - 1) / 180.0))
        z = H[:, j]; t = T[:, j]; p = P[:, j] if P is not None else None
        for (nm, lo, hi, c) in bands:                   # mixed-phase band of each species, where this column has it
            inb = (t >= lo) & (t <= hi)
            if inb.any():
                ax.axhspan(z[inb].min(), z[inb].max(), color=c, alpha=0.07, lw=0)
                ax.text(xmax * 0.6, z[inb].min() + (0.3 if nm == "H2S" else 0.5) * (z[inb].max() - z[inb].min()), nm + " mixed phase",
                        color=c, fontsize=FS_BAND, ha="right", va="center", alpha=0.9)
        if ATREYA:                                       # Atreya ECCM: shaded by pressure, dashed by temperature
            lp = np.log10(np.where(p > 0, p, np.nan)) if p is not None else None
            k_tp = int(np.nanargmin(t)) if np.isfinite(t).any() else len(t) - 1     # troposphere: bottom .. T minimum
            tz, zz = t[:k_tp + 1], z[:k_tp + 1]
            for (sp_, ph_, aT, ap, ar) in ATREYA:
                c_ = acol.get(sp_, "grey"); ar = np.asarray(ar, float)
                if lp is not None and np.isfinite(lp).any():
                    o_ = np.argsort(lp); zp = np.interp(np.log10(ap), lp[o_], z[o_], left=np.nan, right=np.nan)
                    ok = np.isfinite(zp)
                    if ok.sum() > 1:
                        ax.fill_betweenx(zp[ok], xmin, np.maximum(ar[ok], xmin), color=c_, alpha=0.13, lw=0)
                        ax.plot(ar[ok], zp[ok], color=c_, lw=0.8, alpha=0.55)
                        if ok[-1]: ax.plot([xmin, ar[-1]], [zp[-1], zp[-1]], color=c_, lw=1.6, alpha=0.6)
                if len(tz) > 1:
                    o_ = np.argsort(tz); zt = np.interp(aT, tz[o_], zz[o_], left=np.nan, right=np.nan)
                    ok = np.isfinite(zt)
                    if ok.sum() > 1:                    # long thin dashes + x at the base: not the model's ice dashes
                        ax.plot(ar[ok], zt[ok], color=c_, lw=1.0, ls=(0, (9, 4)), alpha=0.7)
                        if ok[-1]: ax.plot([ar[-1]], [zt[-1]], marker="x", ms=7, mew=1.5, color=c_, alpha=0.8)
        for (nm, liq, ice, c, f) in sp:
            L = arr(liq); I = arr(ice)
            if L is not None and (L[:, j] > 0).any(): ax.plot(np.where(L[:, j] > 0, L[:, j] * f, np.nan), z, color=c, lw=1.8, label=nm + " liquid cloud")
            if I is not None and (I[:, j] > 0).any(): ax.plot(np.where(I[:, j] > 0, I[:, j] * f, np.nan), z, color=c, lw=1.8, ls="--", label=nm + " ice")
        S = arr("NH4SH")
        if S is not None and (S[:, j] > 0).any():
            s = S[:, j] * 1e-6                            # written x1e9 (ug/m3) -> g/m3
            ax.plot(np.where(s > 0, s, np.nan), z, color="#8c564b", lw=1.8, ls=":", label="NH4SH (solid, chemistry)")
        ax.set_xscale("log"); ax.set_xlim(xmin, xmax); ax.grid(alpha=0.3, which="both")
        ax.set_xlabel("condensate density [g/m3]", fontsize=FS_LABEL); ax.set_title("lat %+d" % lat, fontsize=FS_TITLE)
        ticks = np.linspace(0, ni - 1, 11).round().astype(int)
        ax.set_ylim(z.min(), z.max()); ax.set_yticks(z[ticks]); ax.set_yticklabels(["%.0f" % t[i] for i in ticks])
        if ax is axs[0]: ax.set_ylabel("temperature [K]", fontsize=FS_LABEL)
        if p is not None:
            ax2 = ax.twinx(); ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(z[ticks])
            ax2.set_yticklabels(["%.3g" % p[i] for i in ticks], fontsize=FS_TICK)
            if ax is axs[-1]: ax2.set_ylabel("pressure [bar]", fontsize=FS_LABEL)
        ax.tick_params(labelsize=FS_TICK)
    h, l = [], []
    for ax in axs:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l: h.append(hh); l.append(ll)
    if ATREYA:
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        h.append(Patch(facecolor="grey", alpha=0.25, edgecolor="none")); l.append("Atreya ECCM by pressure (shaded, bar = base)")
        h.append(Line2D([0], [0], color="grey", lw=1.0, ls=(0, (9, 4)), marker="x")); l.append("Atreya ECCM by temperature (long dashes, x = base)")
        fig.text(0.5, 0.005, "Atreya: " + ATREYA_LABEL + "; densities g/L x 1000, upper limits; colour = species", ha="center", va="bottom", fontsize=FS_BAND)
    if h: fig.legend(h, l, loc="lower center", ncol=min(len(l), 5), fontsize=FS_LEGEND, frameon=False, bbox_to_anchor=(0.5, 0.03 if ATREYA else 0.0))
    fig.suptitle(TITLE, fontsize=FS_SUPTITLE)
    fig.subplots_adjust(left=0.065, right=0.93, bottom=0.2, top=0.88, wspace=0.34)
    return python_view.figure_to_image(fig)
'''

ATREYA_DEFAULT = {"Jupiter": "1x", "Saturn": "5x", "Uranus": "SN16_N5", "Neptune": "SN16_N5"}
def atreya_decks(planet, case):
    """the decks of one planet / case from literature/atreya_eccm.csv: [(species, phase, [T], [p], [rho g/m3])], label"""
    import csv
    if case in (None, "", "none"): return [], ""
    case = ATREYA_DEFAULT.get(planet) if case == "auto" else case
    f = os.environ.get("ATREYA_ECCM") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "literature", "atreya_eccm.csv")
    if not case or not os.path.exists(f): print("no Atreya overlay (%s, %s)" % (case, f)); return [], ""
    decks, src = {}, ""
    with open(f) as h:
        for r in csv.DictReader(h):
            if r["planet"] != planet or r["case"] != case: continue
            src = r["source"]; d = decks.setdefault((r["species"], r["phase"]), ([], [], []))
            d[0].append(float(r["T_K"])); d[1].append(float(r["p_bar"])); d[2].append(float(r["rho_g_per_L"]) * 1000.0)
    if not decks: print("no Atreya case %s for %s in %s" % (case, planet, f)); return [], ""
    return [(k[0], k[1], v[0], v[1], v[2]) for k, v in decks.items()], "%s, case %s" % (src, case)

def build(fname, lats, size=(1800, 900), out=None, atreya="auto"):
    from paraview.simple import LegacyVTKReader, CreateView, Show, Render, SaveScreenshot, GetActiveSource, AssignViewToLayout
    base = os.path.basename(fname)
    planet = base.split("_")[0]
    it = base.rsplit("_", 1)[-1].split(".")[0]
    run = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(fname))))
    title = "%s %s, iteration %s: cloud and ice densities (liquid solid, ice dashed, NH4SH dotted; shaded: mixed-phase band)" % (planet, run, it)
    src = LegacyVTKReader(registrationName=base, FileNames=[fname])
    view = CreateView("PythonView")
    ad, al = atreya_decks(planet, atreya)
    view.Script = (VIEW_SCRIPT.replace("@LATS@", repr(list(lats))).replace("@PLANET@", repr(planet)).replace("@TITLE@", repr(title))
                   .replace("@ATREYA@", repr(ad)).replace("@ATREYA_LABEL@", repr(al)))
    view.ViewSize = list(size)
    try: AssignViewToLayout(view)
    except Exception: pass
    Show(src, view); Render(view)
    if out:
        SaveScreenshot(out, view, ImageResolution=list(size)); print("wrote", out)
    return view

def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file"); ap.add_argument("--lat", help="latitudes, e.g. 0,45,75 (default per planet, see above)")
    ap.add_argument("--out"); ap.add_argument("--size", default="1800x900")
    ap.add_argument("--atreya", default="auto", help="Atreya ECCM case to overlay (auto = per planet, none = off)")
    a = ap.parse_args()
    if not a.lat: a.lat = ",".join(str(x) for x in default_lats(a.file))
    out = a.out or os.path.splitext(a.file)[0] + "_cloudprofile.png"
    build(a.file, [float(x) for x in a.lat.split(",")], tuple(int(x) for x in a.size.split("x")), out, a.atreya)

def default_lats(fname):
    return [0, 45, 75]

if __name__ == "__main__":
    if "paraview.simple" in sys.modules and len(sys.argv) <= 1:     # Python Shell > Run Script inside the GUI
        from paraview.simple import GetActiveSource
        s = GetActiveSource(); f = None
        if s is not None and hasattr(s, "FileNames") and s.FileNames and "_zonal_" in s.FileNames[0]: f = s.FileNames[0]
        f = f or os.environ.get("CLOUDPROFILE_FILE") or DEFAULT_FILE
        if not f: raise SystemExit("open a <Planet>_zonal_*.vtk first (or set CLOUDPROFILE_FILE / DEFAULT_FILE)")
        build(f, default_lats(f), atreya=os.environ.get("CLOUDPROFILE_ATREYA", "auto"))
    else:
        main()
