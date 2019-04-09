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

#include "LinearFitThermFunScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermFunSchemes 
{
    defineTypeNameAndDebug(LinearFitThermFunScheme, 0);
    addToRunTimeSelectionTable(thermFunScheme, LinearFitThermFunScheme, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermFunSchemes::LinearFitThermFunScheme::LinearFitThermFunScheme
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermFunScheme(name),
    //thermFunScheme(name + "=thermFun(T)"),
    //thermFunScheme(name, T),
    //thermFunScheme(name + "=thermFun(T)", T),
    c1_(dict.lookup("c1")),
    c2_(dict.lookup("c2")),
    F0_(dict.lookup(name + "0")),
    F_
    (
        IOobject
        (
            name,
            T.time().timeName(),
            T.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        F0_
        //dimensionedScalar(F.name(), F0_.dimensions(), F0_.value())
    )
{
    //Info<< "LinearFit, F.dimensions() = " << F.dimensions() << endl;
    //Info<< "constantFit, F0_.name()       = " << F0_.name() << endl;
    //Info<< "constantFit, F0_.dimensions() = " << F0_.dimensions() << endl;
    //Info<< "constantFit, F0_.value()      = " << F0_.value() << endl;
    //Info<< "LinearFit, F_.name()       = " << F_.name() << endl;
    //Info<< "LinearFit, F_.dimensions() = " << F_.dimensions() << endl;
    //Info<< "LinearFit, F_[0]           = " << F_[0] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermFunSchemes::LinearFitThermFunScheme::calcThermFun
(
    const volScalarField& T,
    const dimensionedScalar& T0
)
{
    //Info<< "executing LinearFitThermFunScheme" << endl;
    F_ = F0_ * (c1_ + c2_ * (T - T0));
}


// ************************************************************************* //
