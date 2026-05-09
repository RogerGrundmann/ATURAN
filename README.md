# ATJUP

ATJUP (Jupiter Atmospheric Model) based on the numerical treatment of ATOM climate model.

## Concise code description:

Expanding the code in a hydrogene and helium atmosphere for the ammonia transport equation.

- Atmosphere Jupiter General Circulation Model ( ATJUP ) applied to laminar flow
- Program for the computation of Jupiter-atmospherical circulating flows in a spherical shell
- Finite-Difference-Scheme for the solution of the 3D Navier-Stokes equations
  with 6 transport equations to describe the water vapour, cloud water, cloud ice and NH3 vapour, NH3 cloud and NH3 ice
- 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
- Poisson equation for the pressure solution in an outer iterational loop
- Temperature distribution given as a parabolic distribution from pole to pole, zonally constant
- Water and NH3 vapour distribution given by Clausius-Claperon equation for the partial pressure
- Water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model
- Two-Category-Ice-Scheme for cold clouds applying parameterization schemes provided by the COSMO code (German Weather Forecast)
- Rain and snow precipitation solved by column equilibrium applying the diagnostic equations
- 2,4 million grid points (360 x 180 x 40)
- Computer time on a laptop approximately 15 min


## Work in progress:

The plot shows a zonal view of the NH3-ice distribution in the background. In the foreground, the horizontal closed white lines surround the water-cloud-ice and the vertical lines show the location of five circulation cells north and south of the equator. The colour of the latter indicate the vertical u-velocity component. Expanding the number of cells is a matter of copying the existing ones according to the measured w-velocity components by the Voyager(1979), Cassini (2000) and HST (Hubble 2015) missions. The measured temperature/pressure distributions by these missions serve as boundary conditions. The vertical extension reaches from
p = 10e6 Pa (T = 350 K) to p = 10e4 Pa (T = 110 K).

![Jupiter NH3-water-ice-clouds in circulation cells](Jupiter_zonal.png)


## Author

Code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz (roger.grundmann@web.de)

