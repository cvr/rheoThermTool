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

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(multiMode, 0);
    addToRunTimeSelectionTable(constitutiveEq, multiMode, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::multiMode::multiMode
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
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    models_()
{
    PtrList<entry> modelEntries(dict.lookup("models"));
    models_.setSize(modelEntries.size());

    forAll (models_, modelI)
    {
        models_.set
        (
            modelI,
            constitutiveEq::New
            (
                word(name + modelEntries[modelI].keyword()),
                U,
                phi,
		T,
                modelEntries[modelI].dict()
            )
        );
    }

    // Density is treated as the average of the density field of the
    // several models (is it really?).
    rhoRef() = models_[0].rho();
    
    // VE Solvent viscosity is treated as the sum of the several etaS fields.
    etaSRef() = models_[0].etaS();
    
    int cnt(1);
    for (label i = 1; i < models_.size(); i++)
    {
        rhoRef() += models_[i].rho();
        etaSRef() += models_[i].etaS();
        cnt++;
    }
    rhoRef() /= cnt;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTau(const volVectorField& U)
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTau(U);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTau(U);
    }

    return divMatrix;
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTauS(const volVectorField& U, const volScalarField& alpha)
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTauS(U, alpha);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTauS(U, alpha);
    }

    return divMatrix;
}


Foam::tmp<Foam::volSymmTensorField> Foam::constitutiveEqs::multiMode::tau() const
{
    tau_ *= 0;

    for (label i = 0; i < models_.size(); i++)
    {
        tau_ += models_[i].tau();
    }

    return tau_;
}

void Foam::constitutiveEqs::multiMode::rhoMulti()
{
    // Density is treated as the average of the density field of the
    // several models (is it really?).
    rhoRef() = models_[0].rho();
    int cnt(1);
    for (label i = 1; i < models_.size(); i++)
    {
        rhoRef() += models_[i].rho();
        cnt++;
    }
    rhoRef() /= cnt;
}

void Foam::constitutiveEqs::multiMode::etaSMulti()
{
    // VE Solvent viscosity is treated as the sum of the several etaS fields.
    etaSRef() = models_[0].etaS();
    for (label i = 1; i < models_.size(); i++)
    {
        etaSRef() += models_[i].etaS();
    }
}

void Foam::constitutiveEqs::multiMode::correct()
{
    forAll (models_, i)
    {
        Info<< "Model mode "  << i+1 << endl;
        models_[i].correct();
    }

    tau();
    rhoMulti();
    etaSMulti();
}

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::multiMode::energyExtraTerms()
{
    volScalarField energyExtraTerms_ = models_[0].energyExtraTerms()();
    for (label i = 1; i < models_.size(); i++)
    {
        energyExtraTerms_ += models_[i].energyExtraTerms()();
    }
    return tmp<volScalarField>(energyExtraTerms_);
}

// ************************************************************************* //
