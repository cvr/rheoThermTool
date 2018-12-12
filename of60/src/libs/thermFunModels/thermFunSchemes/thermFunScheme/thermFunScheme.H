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

Class
    thermFunScheme

Description
    Abstract base class from which each specific thermal function is derived.

\*---------------------------------------------------------------------------*/

#ifndef thermFunScheme_H
#define thermFunScheme_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thermFunScheme Declaration
\*---------------------------------------------------------------------------*/

class thermFunScheme
{
private:

    // Private data

        //- Name
        const word name_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermFunScheme(const thermFunScheme&);

        //- Disallow default bitwise assignment
        void operator=(const thermFunScheme&);
        

public:

    //- Runtime type information
    TypeName("thermFunScheme");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermFunScheme,
            dictionary,
            (
                const word& name,
                const volScalarField& T,
                const dictionary& dict
            ),
            (name, T, dict)
        );


    // Constructors

        //- Construct from components
        thermFunScheme
        (
            const word& name
        )
        :
            name_(name)
        {}

    // Selectors

        //- Return a reference to the selected thermoFunScheme
        static autoPtr<thermFunScheme> New
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );

    // Destructor

        virtual ~thermFunScheme()
        {}

    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }
        
        //- Return true if thermFunScheme is constant
        virtual bool isConstant() const = 0;

        //- Return the F field
        virtual tmp<volScalarField> F() const = 0;
        
        //- Return F0
        virtual const dimensionedScalar F0() const = 0;
        
        //- Apply thermFun to compute F_ 
        virtual void calcThermFun
        (
            const volScalarField&,
            const dimensionedScalar&
        ) = 0;

        //- Apply thermFun to compute F_ and get F_
        virtual tmp<volScalarField> applyThermFun
        (
            const volScalarField&,
            const dimensionedScalar&
        ) = 0;

        //- Return the F field at element faces
        virtual tmp<surfaceScalarField> Fface() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
