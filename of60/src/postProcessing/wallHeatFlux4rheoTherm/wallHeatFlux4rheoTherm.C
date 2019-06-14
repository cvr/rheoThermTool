/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Changelog:
2019-03-27: Carlos Veiga Rodrigues: adaptation to rheoThermTool (OpenFOAM-6)
2018-06-15: Bruno Santos @ FSD blueCAPE Lda: Adapted to OpenFOAM 5.x.
2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.
2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»
2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»
2010-08-02 Eelco van Vliet: 1st public version of wallHeatFluxIncompressible:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wallHeatFlux4rheoTherm

Description
    Computes the heat flux for all patches as the boundary field of a
    volScalarField, printing the integrated flux for all wall patches.
    Based on wallHeatFlux with changes to allow it for rheoThermTool solvers.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        Info<< "Reading field T\n" << endl;
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        Info << "Reading field kt\n" << endl;
        volScalarField kt
        (
            IOobject
            (
                "kt",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        Info << "Reading field cp\n" << endl;
        volScalarField cp
        (
            IOobject
            (
                "cp",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        Info << "Reading field rho\n" << endl;
        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh
        );

        Info << "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info << "Creating field alphat\n" << endl;
        volScalarField alpha
        (
            IOobject
            (
                "alpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("alpha", dimLength*dimLength/dimTime, 0)
        );
        alpha = kt / rho / cp;

        #include "createPhi.H"

        //label pRefCell = 0;
        //scalar pRefValue = 0.0;
        //setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

        //Correct boundaries only after all fields are loaded in
        T.correctBoundaryConditions();

        // // Thermal expansion coefficient [1/K]
        // dimensionedScalar Pr(laminarTransport.lookup("Pr"));
        // // Turbulent Prandtl number
        // dimensionedScalar Prt(laminarTransport.lookup("Prt"));
        // // Heat capacity
        // dimensionedScalar Cp0(laminarTransport.lookup("Cp0"));
        // // Fluid density
        // dimensionedScalar rho0(laminarTransport.lookup("rho0"));


        /*
        if
        (
            !IOobject("alphat", runTime.timeName(), mesh)
                .typeHeaderOk<volScalarField>(true)
        )
        {
            Info<< "\nCalculating thermal diffusivity " << endl;
            alphat = turbulence->nut()/Prt;
            alphat.correctBoundaryConditions();
        }
        else
        {
            Info<< "\nRead thermal diffusivity alpha" << endl;
        }
        */
        
        surfaceScalarField snGradT
        (
            IOobject
            (
                "snGradT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("snGradT", dimTemperature/dimLength, scalar(0))
        );

        snGradT = fvc::snGrad(T);

        surfaceScalarField heatFlux = fvc::interpolate(kt) * snGradT;

        const surfaceScalarField::Boundary& patchGradT =
                 snGradT.boundaryField();

        const surfaceScalarField::Boundary& patchHeatFlux =
                 heatFlux.boundaryField();

        Info<< "\nWall heat fluxes " << endl;
        forAll(patchHeatFlux, patchi)
        {
           if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
            {
                Info<< mesh.boundary()[patchi].name()
                    << ": Total "
                    << sum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << " [W] over "
                    << sum
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [m2] ("
                    << sum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )/
                       sum 
                       (
                           mesh.magSf().boundaryField()[patchi]
                       )
                    << " [W/m2])"
                    << endl;
            }
      }
      Info<< endl;

      volScalarField wallHeatFlux
        (
            IOobject
            (
                "wallHeatFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
        );

      volScalarField wallGradT
        (
            IOobject
            (
                "wallGradT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("wallGradT", snGradT.dimensions(), 0.0)
        );

      forAll(wallHeatFlux.boundaryField(), patchi)
      {
         wallHeatFlux.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];
      }

      forAll(wallGradT.boundaryField(), patchi)
      {
         wallGradT.boundaryFieldRef()[patchi] = patchGradT[patchi];
      }

      volVectorField gradT = fvc::grad(T);

      Info<< "Writing field '" << wallGradT.name() << "' to " << runTime.timeName() << "/" <<endl;
      wallGradT.write();
      Info<< "Writing field '" << gradT.name() << "' to " << runTime.timeName() << "/" <<endl;
      gradT.write();
      Info<< "Writing field '" << wallHeatFlux.name() << "' to " << runTime.timeName() << "/" <<endl;
      wallHeatFlux.write();
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
