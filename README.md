# rheoThermTool

A toolbox to simulate of thermo-rheological fluid flows in OpenFOAM®.
It contains models for Viscoelastic fluids (VE) and Generalized Newtonian Fluids (GNF)
and the ability to change fluid properties according to temperature-dependent functions.

The `rheoThermTool` software has a framework which allows a user to choose if a quantity
(e.g. reference viscosity or thermal conductivity) is subjected to a thermal law
such as the Arrhenius equation.

The initial implementation for the non-newtonian fluid models
was based from the [rheoTool](https://github.com/fppimenta/rheoTool) toolbox.
The thermo-rheological aspects of the fluid models follows
the [Peters & Baaijens, 1997, J. Non-Newtonian Fluid Mech.](http://doi.org/10.1016/S0377-0257(96)01511-X)
framework for both the constitutive and energy equations.

This software was developed and tested using the [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6) open source CFD toolbox from OpenFOAM Foundation.


## Libraries compiled by `rheoThermTool`

The default location for the libraries of `rheoThermTool`
is the `$FOAM_USER_LIBBIN` directory.  The libraries are:

* `libThermFunModels.so`
* `libRheoThermConstitutiveEq.so`
* `libBCRheoThermTool.so`
* `libpostProcessingRheoThermTool.so`
* `libgaussDefCmpwConvectionSchemes.so`

To change the default location for the libraries,
at compilation time one may export the `FOAM_USER_LIBBIN`
environment variable and point it to a different location.
For example, if one requires the libraries are installed
at folder `$FOAM_LIBBIN` instead:
```sh
export WM_PROJECT_USER_DIR=$WM_PROJECT_DIR
export FOAM_USER_LIBBIN=$FOAM_LIBBIN

cd of60/src
./Allwmake 2>&1 | tee -a logAllwmake.out
```

## Change default location for `rheoThermTool` libraries and applications

To change the default location of libraries and applications,
at compilation time one may export the `WM_PROJECT_USER_DIR`,
`FOAM_USER_LIBBIN` and `FOAM_USER_APPBIN` environment variables,
pointing these to a different locations.
For example, to have `rheoThermTool` in a folder named
`/opt/rheoThermTool`:
```sh
mkdir -p /opt/rheoThermTool/lib
mkdir -p /opt/rheoThermTool/bin

export WM_PROJECT_USER_DIR=/opt/rheoThermTool
export FOAM_USER_LIBBIN=/opt/rheoThermTool/lib
export FOAM_USER_APPBIN=/opt/rheoThermTool/bin

cd of60/src
./Allwmake 2>&1 | tee -a logAllwmake.out
```













## Copyright notice

Copyright (C) 2018 by Carlos Veiga Rodrigues <carlos.rodrigues@fe.up.pt>
and Alexandre Afonso Afonso <aafonso@fe.up.pt>. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For more details consult the GNU General Public License at:
<http://www.gnu.org/licenses/gpl.html>.

