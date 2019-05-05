/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Changelog:
2019-05-05: Carlos Veiga Rodrigues: first version (OpenFOAM-6)

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
    elasticWorkDissipationRatio4RheoTherm

Description
    Computes the ratio between elastic work and dissipation.
    This software is a part of the rheoThermTool Solver.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "constitutiveEq.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "IOdictionary.H"
#include "timeSelector.H"

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
        
        //if (runTime.timeName() == "0")
        //{
        //    Info<< "skipping ...\n" <<endl;
        //    continue;
        //}
        
        mesh.readUpdate();

        IOobject isUPresent
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        );
        if (!isUPresent.typeHeaderOk<volVectorField>(true))
        {
            Info<< "No U field, skipping ...\n" <<endl;
            continue;
        }
        
        IOobject isPhiPresent
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );
        if (!isPhiPresent.typeHeaderOk<surfaceScalarField>(true))
        {
            Info<< "No phi field, skipping ...\n" <<endl;
            continue;
        }
        
        IOobject isTauPresent
        (
            "tau",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );
        if (!isTauPresent.typeHeaderOk<volSymmTensorField>(true))
        {
            Info<< "No tau field, skipping ...\n" <<endl;
            continue;
        }
        
        IOobject isEtaPPresent
        (
            "etaP",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );
        if (!isEtaPPresent.typeHeaderOk<volScalarField>(true))
        {
            Info<< "No etaP field, skipping ...\n" <<endl;
            continue;
        }
        
        IOobject isLambdaPresent
        (
            "lambda",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );
        if (!isLambdaPresent.typeHeaderOk<volScalarField>(true))
        {
            Info<< "No lambda field, skipping ...\n" <<endl;
            continue;
        }

        /*
        Info<< "Reading field p\n" << endl;
        volScalarField p
        (   
            IOobject
            (   
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
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
        */

        Info << "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        //#include "createPhi.H"
        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
            //fvc::flux(U)
        );

        //Info<< "Reading field T\n" << endl;
        //volScalarField T
        //(
        //    IOobject
        //    (
        //        "T",
        //        runTime.timeName(),
        //        mesh,
        //        IOobject::MUST_READ,
        //        IOobject::NO_WRITE
        //    ),
        //    mesh
        //);

        Info<< "Preparing IO for constant/constitutiveProperties\n" << endl;
        IOdictionary cttProperties
	(
	    IOobject
	    (
		"constitutiveProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE,
		false
	    )
	);

        word typeName (cttProperties.subDict("parameters").lookup("type"));
	Info<< "Constitutive equation type: " << typeName << endl << endl;

        //dimensionedScalar elastEnergDiss
        //(
        //    cttProperties.subDict("parameters").lookup("elastEnergDiss")
        //);
	
        //dimensionedScalar T0 = cttProperties.subDict("thermophysicalProperties").lookup("T0");

        /*
        Info<< "\nCreating model(s) for constitutive equation" << endl;
        constitutiveModel constEq(U, phi, T);
	Info<< "Constitutive equation type: " << constEq.eqTypeName() << endl << endl;
        */

        Info << "Reading field tau\n" << endl;
	volSymmTensorField tau
	(
	    IOobject
	    (
		"tau",
                runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        //Info << "Reading field etaS\n" << endl;
	//volScalarField etaS
	//(
	//    IOobject
	//    (
	//	"etaS",
	//	runTime.timeName(),
	//	mesh,
	//	IOobject::MUST_READ,
	//	IOobject::NO_WRITE
	//    ),
	//    mesh
	//);

        Info << "Reading field etaP\n" << endl;
	volScalarField etaP
	(
	    IOobject
	    (
		"etaP",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        Info << "Reading field lambda\n" << endl;
	volScalarField lambda
	(
	    IOobject
	    (
		"lambda",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    ),
	    mesh
	);

        //Info << "Reading field kt\n" << endl;
        //volScalarField kt
        //(
        //    IOobject
        //    (
        //        "kt",
        //        runTime.timeName(),
        //        mesh,
        //        IOobject::READ_IF_PRESENT,
        //        IOobject::NO_WRITE
        //    ),
        //    mesh
        //);
        
        //Info << "Reading field cp\n" << endl;
        //volScalarField cp
        //(
        //    IOobject
        //    (
        //        "cp",
        //        runTime.timeName(),
        //        mesh,
        //        IOobject::READ_IF_PRESENT,
        //        IOobject::NO_WRITE
        //    ),
        //    mesh
        //);
        
        //Correct boundaries only after all fields are loaded in
        //T.correctBoundaryConditions();

        bool validConstitutiveModel = true;
        
        volScalarField eDWr
        (
            IOobject
            (
                "eDWr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        
        volScalarField eDWr_den = (tau && fvc::grad(U)) + etaP/lambda*fvc::div(phi);
        /*
        forAll(eDWr_den.boundaryField(), patchi)
        {
            forAll (eDWr_den.boundaryField()[patchi], facei)
            {
                if ( eDWr_den.boundaryField()[patchi][facei] > -ROOTSMALL
                  && eDWr_den.boundaryField()[patchi][facei] < ROOTSMALL )
                {
                    if (eDWr_den.boundaryField()[patchi][facei] >= 0)
                    {
                        eDWr_den.boundaryFieldRef()[patchi][facei] = ROOTSMALL;
                    }
                    else
                    {
                        eDWr_den.boundaryFieldRef()[patchi][facei] = -ROOTSMALL;
                    }
                }
            }
        }
        */
        forAll(eDWr_den, celli)
        {
            if (eDWr_den[celli] > -ROOTSMALL && eDWr_den[celli] < ROOTSMALL)
            {
                if (eDWr_den[celli] >= 0)
                {
                    eDWr_den[celli] = ROOTSMALL;
                }
                else
                {
                    eDWr_den[celli] = -ROOTSMALL;
                }
            }
        }
        /*
        forAll(lambda.boundaryField(), patchi)
        {
            forAll (lambda.boundaryField()[patchi], facei)
            {
                if ( lambda.boundaryField()[patchi][facei] > -ROOTSMALL
                  && lambda.boundaryField()[patchi][facei] < ROOTSMALL )
                {
                    if (lambda.boundaryField()[patchi][facei] >= 0)
                    {
                        lambda.boundaryFieldRef()[patchi][facei] = ROOTSMALL;
                    }
                    else
                    {
                        lambda.boundaryFieldRef()[patchi][facei] = -ROOTSMALL;
                    }
                }
            }
        }
        */
        forAll(lambda, celli)
        {
            if (lambda[celli] > -ROOTSMALL && lambda[celli] < ROOTSMALL)
            {
                if (lambda[celli] >= 0)
                {
                    lambda[celli] = ROOTSMALL;
                }
                else
                {
                    lambda[celli] = -ROOTSMALL;
                }
            }
        }

        if (word("PTTexp") == typeName)
        {
            dimensionedScalar epsilon
            (
                cttProperties.subDict("parameters").lookup("epsilon")
            );
        
            eDWr = Foam::exp(epsilon*lambda/etaP*tr(tau)) * tr(tau)
                / (2.0 * lambda * eDWr_den);
        }
        else if (word("PTTlinear") == typeName)
        {
            dimensionedScalar epsilon
            (
                cttProperties.subDict("parameters").lookup("epsilon")
            );
            
            eDWr = (epsilon*lambda/etaP*tr(tau) + 1.0) * tr(tau)
                / (2.0 * lambda * eDWr_den);
        }
        else if (word("Giesekus") == typeName)
        {
            dimensionedScalar alpha
            (
                cttProperties.subDict("parameters").lookup("alpha")
            );
            
            eDWr = (alpha * lambda/etaP * (tau && tau) + tr(tau))
                / (2.0 * lambda * eDWr_den);
        }
        else if (word("Oldroyd-B") == typeName)
        {
            eDWr = tr(tau) / (2.0 * lambda * eDWr_den);
        }
        else
        {
            validConstitutiveModel = false;
            Info<< "Constitutive Equation Model not one of:"
                << "\n\t'Giesekus'"
                << "\n\t'Oldroyd-B'"
                << "\n\t'PTTexp'"
                << "\n\t'PTTlinear'"
                //<< "\n... aborting"
                <<endl;
            //return 1;
        }

        forAll(eDWr.boundaryField(), patchi)
        {
            forAll (eDWr.boundaryField()[patchi], facei)
            {
                eDWr.boundaryFieldRef()[patchi][facei] = 0;
            }
        }

        if (validConstitutiveModel)
        {
            Info<< "Writing field '"
                << eDWr.name()
                << "' to "
                << runTime.timeName()
                << "/" <<endl;
            eDWr.write();
        }
        else
        {
            Info<< "Skipping for time " << runTime.timeName() <<endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
