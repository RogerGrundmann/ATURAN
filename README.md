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

---

## Usage

### Command-line

```bash
./cli/uran uranus/config_aturan.xml
```

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

## Author

Roger Grundmann — roger.grundmann@web.de
