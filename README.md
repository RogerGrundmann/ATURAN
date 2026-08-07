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
environment variable to enable. All of them fill diagnostic arrays. Only one of them can feed back
into the temperature equation — `ATURAN_RAD_COUPLING`, added after ATJUP's and ATSAT's; precipitation
and turbulence remain diagnostic-only here.

**Radiation** — grey two-stream, shared `Radiation.h`

| Variable | Default | Effect |
|----------|---------|--------|
| `ATURAN_RADIATION` | 0 | run the solve; fills `Q_rad`, `radiation`, `epsilon` |
| `ATURAN_RAD_COUPLING` | 0.0 | add `Q_rad` to `rhs_t` as `Q_rad·L_rad/(ρ·cp·u_0·t_ref)`. **1.0 is the physically correct value, not a starting point** — see *Known limitations* for why it is invisible at that setting and what the sweep measured |
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
| `ATURAN_PRESS_SOLVER` | 0 | 0 = this model's own serial Gauss-Seidel `computePressure()`; 1 = the shared red-black `PressureSolver<Planet>` |
| `ATURAN_STEADY` | 1 | steady-state query in the report |
| `ATURAN_LOCAL_RHO`, `ATURAN_METRIC_RADIUS`, `ATURAN_COSTHE_ABS`, `ATURAN_PDYN_UNITS` | — | legacy/behaviour switches |

**Buoyancy and the hydrostatic split** — all default off; the model is bit-identical with them unset.

| Variable | Default | Effect |
|----------|---------|--------|
| `ATURAN_BUOY_SCALE` | 1.0 | multiplier on the buoyancy in `rhs_u` |
| `ATURAN_BUOY_REF` | 0 | ATSAT's device: subtract the area-weighted horizontal mean in place |
| `ATURAN_HYDRO_SPLIT` | 0 | ATJUP's: carry the buoyancy in `p_hydro`, drop the radial term from `rhs_u`, let the horizontal gradient enter `rhs_v`/`rhs_w` |
| `ATURAN_HYDRO_REF` | 0 | integrate downward from the top instead of up from the deep boundary. Not the intended setting |
| `ATURAN_HYDRO_ND_KM` | 0 | restore the dropped km→m factor in `p_hydro`'s nondimensionalisation, i.e. the pre-`9da831a` behaviour, for attribution only |

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

**That convention is load-bearing, not tidiness: THIS MODEL IS NOT REPRODUCIBLE ABOVE ONE THREAD.**
Two runs of the same binary, same config, both at `OMP_NUM_THREADS=16`, differ — 8 of 13 output
files identical. Same-config run-to-run divergence at a fixed thread count is a race, not rounding.
Measured further:

- 1 thread vs 16 threads: 0 of 13 files match; after 4 iterations 28 of 64 field arrays in the
  zonal slice differ, starting from `u`, `v`, `w`, `t` at **iteration 2** — iteration 1 is clean.
- With every `#pragma omp` disabled: 13 of 13 identical run-to-run **and** identical to the
  1-thread reference. So the cause is OpenMP, not uninitialised memory.
- The site is NOT identified. Serialising files one group at a time makes every subset
  reproducible — the integrator alone, the physics group, the init routines, `BoundaryConditions.h`
  alone — while the full build still diverges. A subset test reduces the parallel work and so can
  hide a race by changing timing; those passes are weak evidence and should not be read as
  clearing those files. **`-fsanitize=thread` is the tool for this, and has not been run yet.**
- `RungeKutta_Uran_Turb.cpp` states its acceptance test as "reproducibility at any thread count".
  The RK stages pass it in isolation; the program does not.

Threading is also a poor trade here. Per time step, 1/2/4/8/16 threads give 8.562/5.876/4.467/
4.320/4.220 s — saturating at **2×**, implying ~46 % of a step is serial. Running N jobs
concurrently at one thread each beats running them sequentially at 16 threads for any N ≥ 3.

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

1. **The photosphere is 77 K too warm, and it is NOT a heating excess — it is unopposed vertical
   redistribution.** Run to 224 iterations with the radiation diagnostic on, the τ=1 photosphere
   settles at **136.45 K against a T_eff(in) of 59.04 K**, emitting **30× the planet's energy
   budget** and still climbing at +0.15 K/iteration. But the column does not gain heat:

   | checkpoint | T(i=0) deep | T(i=20) mid | T(i=40) top | column mean |
   |---|---|---|---|---|
   | 1 | 405.11 | 253.58 | 77.15 | 245.99 |
   | 28 | 346.63 | 246.31 | 139.07 | 243.87 |

   The deep loses 58.5 K, the top gains 61.9 K, and the **mean moves −0.9 %**. Thermal diffusion
   flattens the initial adiabat and nothing anchors the top of the column to the planet's energy
   budget, so the photosphere drifts to roughly the column mean. Until this is fixed, **the opacity
   constants cannot be judged against this model at all**: a photosphere 77 K too warm says nothing
   about kappa. This remains the highest-value open item.

   What has been ruled out by measurement, so it is not re-derived:
   - **Latent heat is not the source.** `Q_Latent`/`Q_Sensible` reach no RHS, and the path that does
     reach `t` — the saturation adjustment — runs as a *sink* here, with ice sublimating throughout
     (ch4_ice 111 → 87, h2o_ice 59 → 53), consistent with the −0.9 % drift.
   - **The buoyancy anomaly and the hydrostatic split do not fix it.** With the km→m factor restored
     *and* the metric radius wired — the configuration both were built to reach — `ATURAN_HYDRO_SPLIT`
     moves OLR/in by 0.34 % alone and 0.007 % on top of the metric. A term that redistributes
     buoyancy *horizontally* is the wrong instrument for a fault that is vertical.
   - **Radiative coupling works, and is far too slow.** See item 3.

2. **That number got worse when a real bug was fixed, and the previous one was not better.** Before
   the methane-viscosity correction this model read 7.591× — an artefact of two errors partly
   cancelling. `mue_ch4` held methane's viscosity in *centipoise* as if it were Pa·s, so the
   mass-weighted `mue_mix` came out ~150× too large, and on this model `mue_mix` sets the species
   diffusivities and hence the diffusive-enthalpy sink in `rhs_t`. That sink, ~150× overweighted, was
   resisting the flattening in item 1. Correcting it unmasked the drift rather than causing it.

3. **`ATURAN_RAD_COUPLING` moves the photosphere the right way and cannot move it far enough.** The
   term is the anchor item 1 is missing, and the sub-cap sweep at nm=224 is monotonic:

   | coupling | T(τ=1) | OLR/in | T(i=40) top |
   |---|---|---|---|
   | 0 (off) | 136.45 | 30.093 | 139.07 |
   | 1.0 | 136.45 | 30.093 | 139.07 |
   | 1e3 | 136.41 | 30.064 | 139.04 |
   | 1e4 | 136.06 | 29.806 | 138.77 |
   | 3e4 | 135.69 | 29.536 | 138.47 |
   | 1e5 | 136.33 | 30.177 | 139.24 |

   **At 1.0 — the physically correct value — the term is live but invisible**: it changes all 92
   output files, and the largest temperature change anywhere is 1.0e-4 K, the output format's own
   resolution. That was predicted before the run from `Q_rad ~ 1e-4 W/m³` and matches, which is what
   certifies the scaling; a *visible* result at 1.0 would have meant a units error. ATSAT's note puts
   the same point at ~1e9 iterations to equilibrate.

   **The 1e5 row is the limiter, not the term.** Its raw tendency is ~1.1 against ATJUP's
   `rad_t_max = 0.5`, so it redistributes by where the cap bites — which is why it breaks the trend
   and warms the top instead of cooling it. Read the sub-cap rows only. Extrapolating those, closing
   the 80 K gap needs a coupling of order 1e6, which is deep in the capped regime: **the term is
   directionally right and cannot reach the answer within the cap at this run length.**

4. **`ATURAN_THERMAL_MASSFLUX` is a measurement instrument, not a fix.** Setting it to 0 removes the
   sink entirely and the model drifts harder (34× at 224 iterations), so the term is load-bearing
   even though its form is questionable: it is a flux times a temperature *gradient magnitude*, with
   an absolute value on one component only, so it cannot change sign to oppose the flattening.

5. **The integrator's temperature floor has been reached in real runs.** `t_min_K()` = 7.5 K was
   once recorded as never engaging, on the evidence of a 2-iteration run; over 224 iterations with
   the pre-fix viscosity it engaged in 20 of 28 checkpoints. It is a guard against a NaN, not a
   dormant one, and a run that touches it is reporting a collapse, not a temperature.

6. **The precipitation and turbulence modules are diagnostic-only here.** They fill their arrays and
   nothing reads them back: `S_precip_*` reaches no RHS and `ATURAN_TURB_COUPLING` defaults to 0.
   Switching either on changes plots, not physics. Radiation is no longer in this list — see item 3.

7. **The model is not reproducible above one thread.** A race, unlocated; see *Usage*. Every
   measurement quoted in this file was taken at `OMP_NUM_THREADS=1`, which is the only setting under
   which the byte-comparisons those measurements rest on are meaningful.

8. **The corrected metric radius is off by default, and the warning attached to it is wrong.**
   `metricRadius()` is the identity unless `ATURAN_METRIC_RADIUS` is set; `rad.z` runs 1..2, so the
   metric puts the surface `L_atm` from the centre instead of R and every horizontal derivative is
   25362/360 ≈ 70× too large. `1c4da64` predicted that correcting it would "quite possibly
   destabilise" the model. Measured at nm=224 with `ATURAN_METRIC_RADIUS=25362`, it does the
   opposite — the continuity residuum falls 0.308 → 0.043 and the meridional wind 12.57 → 0.24 m/s,
   the latter matching the ~47× the correction applies. OLR/in moves 30.093 → 29.581. Whether to
   flip the default is an open decision, but it should be taken against these numbers rather than
   against the warning.

9. **The grey opacity is Jupiter's calibration, not Uranus's.** `C_cia` and `opac_cal` were tuned so
   that Jupiter's photosphere lands at 0.25–0.35 bar. Nothing has recalibrated them here, and the
   per-planet lever is the runtime knob, not a second copy of the constant.

---

## Author

Roger Grundmann — roger.grundmann@web.de
