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
the [Peters & Baaijens (1997, J. Non-Newtonian Fluid Mech.)](http://doi.org/10.1016/S0377-0257(96)01511-X)
framework for both the constitutive and energy equations.

This software was developed and tested using the [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6) open source CFD toolbox from OpenFOAM Foundation.


## Features

Fluid models available:
* VE fluid models: Oldroyd-B, Giesekus, PTT exponential or linear.
* GNF fluid models: Power Law, Carreau-Yasuda and the standard Newtonian stress tensor.

Thermal functions available:
* Williams-Landel-Ferry equation.
* Arrhenius equation.
* Linear fit.
* Constant (i.e. no temperature dependency).

Other Characteristics
* Fully-parallelized solver.
* Generic grids for 2D/3D problems.
* Solvers support either moving or static meshes.
* Transient flow solvers are highly stable regarding pressure-stress-velocity coupling.
* Two stabilization methods for the Momentum Equations based on *Both-Side Diffusion*.
* Log-conformation method to stabilize the Constitutive Equation, available for all VE fluid models.


## To be implemented

* Solver for two-phase flows, using any GNF or VE fluid models for each phase.
* Wider range of tutorials to illustrate the application of the solvers to different problems.


## Compatibility

The development and testing of `rheoThermTool` was mainly
performed in [OpenFOAM v6](https://github.com/OpenFOAM/OpenFOAM-6) from the
[The OpenFOAM Foundation Ltd](https://openfoam.org/). It is, however, expected to
work with both [OpenFOAM v5.x](https://github.com/OpenFOAM/OpenFOAM-5.x)
and [OpenFOAM v4.x](https://github.com/OpenFOAM/OpenFOAM-4.x).
It is also expected to work with the ESI [OpenFOAM®](https://www.openfoam.com)
such as [OpenFOAM® v1806](https://www.openfoam.com/releases/openfoam-v1806).

It was tested on a system running Debian Linux OS, version 9, in an x86_64 architecture and
with the GCC compiler version 6.3.0.


## Third-Party software

The following third-party software is used in `rheoThermTool`:

* [Eigen](http://eigen.tuxfamily.org/).


## Installation

### Requirements

* Compatible and functional version of [OpenFOAM-6](https://github.com/OpenFOAM/OpenFOAM-6) or [OpenFOAM®-v1806](https://develop.openfoam.com/Development/OpenFOAM-plus), already installed.
* Internet connection to download the [Eigen](http://eigen.tuxfamily.org/) C++ linear algebra library.


### Basic installation instructions for OpenFOAM-6

These instructions are appropriate for installation under the user's home path.
For the sake of simplicity, it is assumed that OpenFOAM-6 was previously installed
at `~/OpenFOAM/OpenFOAM-6`. Also, string `user` will be used as a placeholder for the user's name.

1. Define environmental variable `OFDIR` pointing to the OpenFOAM installation:
```sh
user:~ $  OFDIR=~/OpenFOAM/OpenFOAM-6
```

2. Source the environment OpenFOAM variables:
```sh
user:~ $  source $OFDIR/etc/bashrc
```

3. Certify the environment variables were correctly loaded, e.g.:
```sh
user:~ $  echo $WM_PROJECT_INST_DIR
/home/user/OpenFOAM/OpenFOAM-6

user:~ $  echo $WM_PROJECT_USER_DIR
/home/user/OpenFOAM/user-6

user:~ $  echo $FOAM_USER_APPBIN
/home/user/OpenFOAM/user-6/platforms/linux64GccDPInt32Opt/bin

user:~ $  echo $FOAM_USER_LIBBIN
/home/user/OpenFOAM/user-6/platforms/linux64GccDPInt32Opt/lib
```

4. Choose a directory where `rheoThermTool` will be installed.
A good candidate is `$WM_PROJECT_USER_DIR`.
Assuming that `rheoThermTool` will be installed at `$WM_PROJECT_USER_DIR`,
go to that folder and clone the `rheoThermTool` repository:
```sh
user:~ $  cd $WM_PROJECT_USER_DIR

user:~/OpenFOAM/user-6 $  git clone https://github.com/cvr/rheoThermTool

user:~/OpenFOAM/user-6 $  cd rheoThermTool
```
With these instructions, `rheoThermTool` source code and tutorials will be located in `$WM_PROJECT_USER_DIR/rheoThermTool`.

5. Install Eigen and define the `$EIGEN_RHEO` variable with the corresponding path for its libraries:
```sh
user:~/OpenFOAM/user-6/rheoThermTool $  ./downloadEigen

user:~/OpenFOAM/user-6/rheoThermTool $  export EIGEN_RHEO=$(find $WM_PROJECT_USER_DIR/ThirdParty/ -maxdepth 1 -type d -iname "eigen*" | head -n 1)

user:~/OpenFOAM/user-6/rheoThermTool $  echo $EIGEN_RHEO
/home/user/OpenFOAM/user-6/ThirdParty/Eigen3.2.9
```

6. Compile corresponding `rheoThermTool` source code:
```sh
user:~/OpenFOAM/user-6/rheoThermTool $  cd of60/src

user:~/OpenFOAM/user-6/rheoThermTool/of60/src $  ./Allwmake 2>&1 | tee -a log.Allwmake
```

7. Consult the `log.Allwmake` file to understand if there were any installation problems.
```sh
user:~/OpenFOAM/user-6/rheoThermTool/of60/src $  grep -i error log.Allwmake

user:~/OpenFOAM/user-6/rheoThermTool/of60/src $  grep "bin/rheoThermFoam" log.Allwmake
     -lm -o /opt/OpenFOAM/ThirdParty-v1806/rheoThermTool/platforms/linux64GccDPInt32Opt/bin/rheoThermFoam
     -lm -o /opt/OpenFOAM/ThirdParty-v1806/rheoThermTool/platforms/linux64GccDPInt32Opt/bin/rheoThermFoamOptimized
```


### Basic installation instructions for OpenFOAM®-v1806

The instructions are the same as for OpenFOAM-6, with the following differences:

* Environment variable `OFDIR` should point to OpenFOAM®-v1806 installation path,
e.g. `ODFIR=~/OpenFOAM/OpenFOAM-v1806`.
* The corresponding `rheoThermTool` source code is under the `of+1806/src` folder.


### Test installation

Go to your run directory, copy the cylinder tutorial case and run it.
```sh
cd $FOAM_RUN
cp -r $WM_PROJECT_USER_DIR/rheoThermTool/of60/tutorials/cylinder ./

cd cylinder
blockMesh
mirrorMesh -noFunctionObjects -overwrite
rheoThermFoam
paraFoam
```

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

Copyright (C) 2018-2019 by Carlos Veiga Rodrigues <carlos.rodrigues@fe.up.pt>
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

