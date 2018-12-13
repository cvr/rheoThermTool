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

#include "nullThermFunScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermFunSchemes 
{
    defineTypeNameAndDebug(nullThermFunScheme, 0);
    addToRunTimeSelectionTable(thermFunScheme, nullThermFunScheme, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermFunSchemes::nullThermFunScheme::nullThermFunScheme
(
    const word& name,
    const volScalarField& T,
    const dictionary&
)
:
    thermFunScheme(name),
    F0_(dimensionedScalar(word::null, dimless, 0)),
    F_(volScalarField::null())
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermFunSchemes::nullThermFunScheme::calcThermFun
(
    const volScalarField& T,
    const dimensionedScalar& T0
)
{
    Info<< "This is a 'null' thermFunScheme... nothing to do.";
    //F_ = F0_;
}


// ************************************************************************* //
