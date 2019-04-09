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

#include "thermFunModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermFunModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermFunModel::thermFunModel
(
    const word& name,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "constitutiveProperties",
            T.time().constant(),
            T.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    //funPtr_(thermFunScheme::New(word::null, T, subDict("thermophysicalProperties")))
    funPtr_(thermFunScheme::New(name, T, subDict("thermophysicalProperties")))
{
    //Info<< "thermFunModel, F.dimensions() = " << F.dimensions() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermFunModel::calcThermFun
(
    const volScalarField& T,
    const dimensionedScalar& T0
)
{
    funPtr_->calcThermFun(T, T0);
}

tmp<volScalarField> thermFunModel::applyThermFun
(
    const volScalarField& T,
    const dimensionedScalar& T0
)
{
    return funPtr_->applyThermFun(T, T0);
}

tmp<surfaceScalarField> thermFunModel::Fface() const
{
    return funPtr_->Fface();
}

bool thermFunModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
