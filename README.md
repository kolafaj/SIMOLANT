# SIMOLANT

Molecular simulation in 2D. See directory `executables/` for the installation.

![SIMOLANT](sources/screenshot.png)

## Aims

* _In teaching physics and chemistry at elementary and high schools:_ A number of phenomena are shown using a two-dimensional molecular model of matter:

  * Condensation of gas and crystallization of liquid on cooling
  * Melting and evaporation on heating
  * Mixing of fluids and gases
  * Capillary action
  * Crystal defects in motion
  * Gas in a gravitational field
  * Impact of a solid body (crystal) to a wall
  * Vicsek model of a flock of birds

* _In a university course of molecular simulations:_ Basic concepts of statistical thermodynamics and molecular simulations can be elucidated:

  * Ergodic and deterministic dynamic systems
  * Nucleation, Ostwald ripening
  * Molecular dynamics at constant energy / temperature / pressure
  * Monte Carlo at constant energy / temperature / pressure
  * Convergence profiles of quantities
  * Radial distribution function
  * Radial density profile, z-density profile
  * Walls and periodic boundary conditions
  * Flying icecube artifact
  * Expert: export of quantities and statistics, keyboard input

* _Student work_

  * Isotherms, phase diagram
  * [Verification of the Clausius-Clapeyron equation](http://old.vscht.cz/fch/en/tools/kolafa/tul/simenw1.pdf)
  * [Pressure outside a droplet/inside cavity (Kelvin equation)](http://old.vscht.cz/fch/en/tools/kolafa/simenw3.pdf)
  * Surface / interfacial tension, contact angle
  * Second virial coefficient and the equation of state
  * Mean square displacement
  * Diffusivity, Arrhenius plot, activation energy
  * Soft penetrable disks
  * Coalescence of droplets
  * Double-minimum potential which might give the Penrose quasicrystal
  * Phase transition in the Vicsek model

## Original pages
  * https://old.vscht.cz/fch/software/simolant/index-en.html
  * https://old.vscht.cz/fch/software/simolant/index-cz.html (in Czech)

## Licences

Copyright (C) 2025 Jiří Kolafa

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html) (GPLv3), as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

This software was developed using the [GNU suite of tools](https://gcc.gnu.org/onlinedocs/libstdc++/manual/license.html) and the [FLTK cross-platform GUI tools](https://www.fltk.org/doc-1.3/license.html).

The following libraries included in the package for Windows are covered by the respective GPLv3-compatible licences:<br />
zlib1.dll, libwinpthread-1.dll, libpng16-16.dll, libjpeg-8.dll, libgcc_s_seh-1.dll.

## Acknowledgements

* Ivo Nezbeda (the first version of SIMOLANT was written in Pascal for a simulation seminar which was held at the Institute of Chemical Technology in about 1985).
* Karel Matas (FLTK advise).
* Numerous students using this software at the University of Chemistry and Technology, Prague and the Technical University Liberec.
