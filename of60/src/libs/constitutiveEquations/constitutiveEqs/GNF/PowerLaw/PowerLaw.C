/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "PowerLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(PowerLaw, 0);
    addToRunTimeSelectionTable(constitutiveEq, PowerLaw, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::PowerLaw::PowerLaw
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& T,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi, T),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
                "0",
                dimensionSet(1, -1, -2, 0, 0, 0, 0),
                pTraits<symmTensor>::zero  
        ),
        extrapolatedCalculatedFvPatchField<symmTensor>::typeName
    ),
    n_(dict.lookup("n")),
    k_
    (
        IOobject
        (
            "k" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("k"))
    ),
    etaMin_
    (
        IOobject
        (
            "etaMin" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("etaMin"))
    ),
    etaMax_
    (
        IOobject
        (
            "etaMax" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("etaMax"))
    )
{
    rhoRef().rename("rho" + name);
    rhoRef().readOpt() = IOobject::READ_IF_PRESENT;
    rhoRef().writeOpt() = IOobject::AUTO_WRITE;
    rhoRef() = volScalarField
    (
        IOobject
        (
            "rho" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("rho"))
    );

    etaRef().rename("eta" + name);
    etaRef().readOpt() = IOobject::READ_IF_PRESENT;
    etaRef().writeOpt() = IOobject::AUTO_WRITE;
    etaRef() = volScalarField
    (
        IOobject
        (
            "eta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zeroEta", dimPressure*dimTime, 0)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::PowerLaw::correct()
{
    etaRef() = max
    (
        etaMin_,
        min
        (
            etaMax_,
            k_*pow 
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                n_.value() - scalar(1.0)
            )
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::PowerLaw::energyExtraTerms()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "energyExtraTerms",
                U().time().timeName(),
                U().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U().mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}

// ************************************************************************* //
