# ATURAN

Uranus atmospheric general circulation model based on the numerical framework of the ATOM
climate model. Solves the 3-D Navier-Stokes equations in a spherical shell extending from
the deep H₂O / NH₃ cloud decks up through the methane condensation layer that gives Uranus
its blue-green appearance.

Uranus is the coldest planetary atmosphere in the Solar System (49 K at the tropopause)
and the only one whose rotation axis lies essentially in the plane of the ecliptic
(97.77° obliquity). One Uranian pole spends 42 Earth-years pointed at the Sun while the
other is in continuous shadow, and yet the observed temperatures at equator and poles
remain remarkably similar — a long-standing puzzle that motivates a 3-D circulation model
with full radiative-convective forcing. Voyager 2's 1986 flyby remains the only in-situ
encounter, so the available constraints are the Voyager temperature/pressure profile,
ground-based and Hubble imagery, and Spitzer/Herschel infrared spectroscopy.

For all relevant data concerning Uranus and its atmosphere the book *Planetary Sciences*
by Imke de Pater and Jack J. Lissauer was indispensable.

---

## Physics & Numerics

- **Domain:** spherical shell, 41 × 181 × 361 grid points (r × θ × φ), ~2.7 million cells
- **Vertical extent:** 360 km atmospheric shell (~9 km per radial step)
- **Dynamics:** finite-difference discretisation of the 3-D Navier-Stokes equations in spherical coordinates
- **Time integration:** 4th-order Runge-Kutta (inner loop) with a Poisson pressure solver (outer loop)
- **Parallelism:** OpenMP shared-memory threading
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile; reference temperature 76.4 K
  - Equatorial value 60.65 K, polar value 49.15 K — reflecting Uranus's anomalous warm-pole
    equilibrium derived from Voyager 2 IRIS retrievals
  - Boussinesq buoyancy approximation
  - Clausius-Clapeyron / Sanchez-Lavega SVP formulation for saturation vapour pressures
  - Mixed-phase (liquid + ice) saturation adjustment with iterative convergence
- **Microphysics:** two-category ice scheme adapted from the COSMO weather-forecast model
- **Boundary conditions:** zonal-wind templates from QuikSCAT and OSCAR observational
  datasets; temperature/pressure profiles from Voyager 2 (1986) and Hubble Space Telescope
- **Planetary constants:** g = 8.69 m/s², Ω = 1.01 × 10⁻⁴ rad/s (≈ 17.24 h sidereal day),
  reference wind speed 100 m/s

---

## Chemical Species

| Species | Phases modelled | Role on Uranus |
|---------|-----------------|----------------|
| CH₄ | vapour · cloud · ice | Dominant condensable; methane ice-haze deck near the tropopause defines the visible blue-green disc |
| H₂O | vapour · cloud water · cloud ice | Deep tropospheric cloud (~50 bar) |
| NH₃ | vapour · cloud · ice | Mid-troposphere reservoir; sequestered as NH₄SH |
| H₂S | vapour | Required for NH₄SH formation; possibly a dominant volatile in Uranus's deep atmosphere |
| NH₄SH | vapour | Heterogeneous reaction NH₃ + H₂S → NH₄SH at ~230 K |

Methane condensation is the defining cloud process at Uranus's low temperatures and is
fully transported through the chemistry pipeline (RHS, RK4, saturation adjustment,
diffusion mass flux).

---

## Repository Layout

```
ATURAN/
├── planet/          # core model (RHS, RK4, thermodynamics, chemistry, I/O)
├── lib/             # array types, config parser, FFT, utilities
├── cli/             # command-line driver (uran)
├── python/          # Cython bindings (pyaturan)
├── uranus/          # run directory (XML config, observational data, output)
│   ├── oscar/       # OSCAR ocean-current dataset (used as zonal-wind template)
│   └── windspeed/   # QuikSCAT surface wind data
├── tinyxml2/        # vendored XML library
├── param.py         # code-generation script (auto-generates parameter files)
└── Makefile
```

---

## Build

**Dependencies:** C++11 compiler with OpenMP support, Python 3, Cython, NumPy.

```bash
# Generate parameter files and build CLI + Python extension
make

# CLI binary only
make uran

# Python extension only
make python

# Clean
make clean
```

The `param.py` script auto-generates several `.inc` / `.pyx` files that parameterise the
model; it runs automatically as part of the build whenever `param.py` itself changes.

### Shared physics headers

Eleven headers in `planet/` are **byte-identical across ATJUP, ATSAT, ATNEPT and ATURAN**:

```
ATPhys.h  BoundaryConditions.h  ConvectiveAdjustment.h  FluxLimiter.h  ParaViewWriter.h
Precipitation.h  PressureSolver.h  Radiation.h  Reporting.h  SaturationAdjustment.h  Turbulence.h
```

There is no submodule and no symlink holding them together — Synology Drive has silently reverted
a working tree once, and a submodule costs friction on every clone. They are plain copies, and
`planet/SHARED.md5` is what makes a divergence loud:

```bash
make check-shared
```

Editing one means: edit it in one repo, copy it to the other three, regenerate its line in **all
four** manifests, and rebuild each.

**`make check-shared` cannot catch everything, and this is the part to read before trusting it.**
It verifies a repo against *its own* manifest, so two repos holding different copies of the same
header both report OK — a state that has already occurred once. The check that does catch it is a
diff between repos, which is why the checksum lines are kept sorted by filename:

```bash
diff ../ATJUP/planet/SHARED.md5 planet/SHARED.md5
```

---

## Optional modules

Every module is **off by default** and a stock run is unaffected by its presence. Set the
environment variable to enable. All of them fill diagnostic arrays; none feeds back into the
temperature equation on this model — there is no `ATURAN_RAD_COUPLING`, unlike ATJUP.

**Radiation** — grey two-stream, shared `Radiation.h`

| Variable | Default | Effect |
|----------|---------|--------|
| `ATURAN_RADIATION` | 0 | run the solve; fills `Q_rad`, `radiation`, `epsilon` |
| `ATURAN_SOLAR` | 1 | absorbed shortwave channel (only acts with `ATURAN_RADIATION`) |
| `ATURAN_SOLAR_STRENGTH` | 1.0 | scale the absorbed insolation |
| `ATURAN_SW_TAU_PER_BAR` | 1.0 | move the shortwave absorption level |
| `ATURAN_CIA_STRENGTH` | 1.0 | scale the H₂/He collision-induced opacity |
| `ATURAN_OPACITY_STRENGTH` | 1.0 | scale the gas-band and cloud opacity |

**Microphysics and turbulence**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATURAN_PRECIP` | 0 | precipitation scheme; fills `P_*`, `Q_precip`, `S_precip_*` (diagnostic only here — nothing reads `S_precip_*`) |
| `ATURAN_TURB` | 0 | run the closure |
| `ATURAN_TURB_MODEL` | *param* | override `turb_model` (`k_epsilon`, `k_omega`, `k_omega_SST`) |
| `ATURAN_TURB_COUPLING` | 0.0 | feed the eddy viscosity into momentum, heat and species diffusion |
| `ATURAN_CONV_ADJ` | 0 | dry convective adjustment |
| `ATURAN_SATADJ` | *see code* | mirrored saturation adjustment |

**Numerics and experiment knobs**

| Variable | Default | Effect |
|----------|---------|--------|
| `ATURAN_THERMAL_MASSFLUX` | 1.0 | scale the diffusive-enthalpy sink in `rhs_t` — see *Known limitations* |
| `ATURAN_SINTHE_MIN` | 0.0 | env floor on sin θ — **not the value in force**: the integrator uses a hardcoded `sinthe_min = 0.4`, so this accessor is not consulted by default |
| `ATURAN_PRESS_SOLVER` | *see code* | pressure-solver selection |
| `ATURAN_STEADY` | 1 | steady-state query in the report |
| `ATURAN_LOCAL_RHO`, `ATURAN_METRIC_RADIUS`, `ATURAN_COSTHE_ABS`, `ATURAN_PDYN_UNITS` | — | legacy/behaviour switches |

---

## Diagnostics

Printed at the `checkpoint` cadence:

- **printMinMax** — max/min with location for every prognostic and diagnostic field, including the
  radiation trio (`radiation`, `Q_rad`, `emissivity`), the turbulence fields and the precipitation
  fluxes. All read zero when their module is off, and are printed regardless: an absent array and a
  zero array look the same, and only one of them means the writer works.
- **Equatorial column profile** — `i`, `p[bar]`, `T[K]`, `eps`, `netRad`, `Q_rad` from the model top
  down to the deep boundary. This is the diagnostic that localises profile faults; a column *mean*
  cannot. It found a defect on ATNEPT where the top of the domain sat at 15 bar, which every
  column-averaged number had reported as a confident 17.2287 bar photosphere.
- **Photosphere line** (with `ATURAN_RADIATION=1`) — mean OLR against the input budget, the τ=1
  level in bar, and the temperature there beside the blackbody flux it implies. The two temperature
  gaps answer different questions: `T(τ=1) − T_eff(OLR)` is the *scheme's* excess, while
  `T(τ=1) − T_eff(in)` is how far the *column* sits from the planet's energy budget.
- **ParaView** — `Radiation`, `Q_rad_mW_m3` and `Emissivity`, plus the turbulence six and the
  precipitation eleven, in all four views (panorama `.vts`; radial, zonal, longitudinal `.vtk`).

---

## Usage

### Command-line — note the TWO arguments

```bash
./cli/uran <path> <config file name>
./cli/uran . config_aturan.xml            # config in the current directory
```

A single combined path fails with `couldn't load config file inside cUranusModel`. This is the
reverse of ATOM's `cli/atm`, and the same convention ATSAT and ATNEPT use.

```bash
OMP_NUM_THREADS=12 ./cli/uran . config_aturan.xml > run.log 2>&1
```

Runs in this repository's measurements are single-threaded (`OMP_NUM_THREADS=1`) so that results
are bit-reproducible and byte-comparisons between builds mean something.

### Python

```python
import sys
sys.path.insert(0, "uranus")
import pyaturan

model = pyaturan.UranusModel()
model.load_config("uranus/config_aturan.xml")
model.run()
```

Output is written as VTK / VTS files for visualisation in ParaView (panorama, sphere,
radial, zonal, and longitudinal cross-sections).

---

## Known limitations

None of these stops a run; all of them affect what a result means.

1. **There is an unopposed heating excess, and it is the highest-value open item.** Run to 224
   iterations with the radiation diagnostic on, the τ=1 photosphere settles at **136.45 K against a
   T_eff(in) of 59.04 K** — 77 K too warm — emitting **30× the planet's energy budget** and still
   climbing at +0.16 K/iteration. Until this is found, **the opacity constants cannot be judged
   against this model at all**: a photosphere 77 K too warm says nothing about kappa.

2. **That number got worse when a real bug was fixed, and the previous one was not better.** Before
   the methane-viscosity correction this model read 7.591× — an artefact of two errors partly
   cancelling. `mue_ch4` held methane's viscosity in *centipoise* as if it were Pa·s, so the
   mass-weighted `mue_mix` came out ~150× too large, and on this model `mue_mix` sets the species
   diffusivities and hence the diffusive-enthalpy sink in `rhs_t`. That sink, ~150× overweighted, was
   holding the column down against the heating excess above. Correcting it unmasked the excess
   rather than causing it.

3. **`ATURAN_THERMAL_MASSFLUX` is a measurement instrument, not a fix.** Setting it to 0 removes the
   sink entirely and the model runs away harder (34× at 224 iterations), so the term is load-bearing
   even though its form is questionable: it is a flux times a temperature *gradient magnitude*, with
   an absolute value on one component only, so it cannot change sign to oppose a runaway.

4. **The integrator's temperature floor has been reached in real runs.** `t_min_K()` = 7.5 K was
   once recorded as never engaging, on the evidence of a 2-iteration run; over 224 iterations with
   the pre-fix viscosity it engaged in 20 of 28 checkpoints. It is a guard against a NaN, not a
   dormant one, and a run that touches it is reporting a collapse, not a temperature.

5. **The radiation, precipitation and turbulence modules are diagnostic-only here.** They fill their
   arrays and nothing reads them back: there is no `ATURAN_RAD_COUPLING`, `S_precip_*` reaches no
   RHS, and `ATURAN_TURB_COUPLING` defaults to 0. Switching a module on changes plots, not physics.

6. **The grey opacity is Jupiter's calibration, not Uranus's.** `C_cia` and `opac_cal` were tuned so
   that Jupiter's photosphere lands at 0.25–0.35 bar. Nothing has recalibrated them here, and the
   per-planet lever is the runtime knob, not a second copy of the constant.

---

## Author

Roger Grundmann — roger.grundmann@web.de
