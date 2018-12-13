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
    thermFunModel

Description
    Base class which can be directly called by applications to use any of 
    the models derived from class thermFunScheme.

    This class is part of rheoThermTool.

Files
    thermFunModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermFunModel_H
#define thermFunModel_H

#include "thermFunScheme.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thermFunModel Declaration
\*---------------------------------------------------------------------------*/

class thermFunModel
:
    public IOdictionary
{
private:

    // Private data

        //- thermFun Eq
        autoPtr<thermFunScheme> funPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermFunModel(const thermFunModel&);

        //- Disallow default bitwise assignment
        void operator=(const thermFunModel&);


public:

    //- Runtime type information
    TypeName("thermFunModel");


    // Constructors

        //- Construct from components
        thermFunModel(const word& name, const volScalarField& T);

    // Destructor

        virtual ~thermFunModel()
        {}


    // Member Functions

        //- Return true if thermFunScheme is active
        virtual bool active() const { return funPtr_->active(); }

        //- Return true if thermFunScheme is constant
        virtual bool isConstant() const { return funPtr_->isConstant(); }
        
        //- Get field F
        virtual tmp<volScalarField> F() const { return funPtr_->F(); }

        //- Get field F0
        virtual const dimensionedScalar F0() const { return funPtr_->F0(); }

        //- Apply thermFun to compute F_ 
        virtual void calcThermFun
        (
            const volScalarField& T,
            const dimensionedScalar& T0
        );
      
        //- Apply thermFun to compute F_ and get F_
        virtual tmp<volScalarField> applyThermFun
        (
            const volScalarField& T,
            const dimensionedScalar& T0
        );

        //- Return the F field at element faces
        virtual tmp<surfaceScalarField> Fface() const;

        //- Read thermFunProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
