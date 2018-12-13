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
    nullThermFunScheme

Description
    Null class for inexistent parameters or for initialization.
    
    This class is part of rheoThermTool.

Files
    nullThermFunScheme.C

\*---------------------------------------------------------------------------*/

#ifndef nullThermFunScheme_H
#define nullThermFunScheme_H

#include "thermFunScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermFunSchemes
{ 

/*---------------------------------------------------------------------------*\
                           Class nullThermFunScheme Declaration
\*---------------------------------------------------------------------------*/

class nullThermFunScheme
:
    public thermFunScheme
{
    // Private data
        
        //- Reference F value (declare before F)
        dimensionedScalar F0_;    

        //- Computed field F (declare after F0)
        volScalarField F_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        nullThermFunScheme(const nullThermFunScheme&);

        //- Disallow default bitwise assignment
        void operator=(const nullThermFunScheme&);


public:

    //- Runtime type information
    TypeName("null");

    // Constructors

        //- Construct from components
        nullThermFunScheme
        (
            const word& name,
            const volScalarField& T,
            const dictionary&
        );


    // Destructor

        virtual ~nullThermFunScheme()
        {}


    // Member Functions
        
        //- Return true if thermFunScheme is active
        virtual bool active() const { return false; }
        
        //- Return true as this thermFun is constant
        virtual bool isConstant() const { return true; }

        //- Return field f
        virtual tmp<volScalarField> F() const { return F_; }
        
        //- Return F0
        virtual const dimensionedScalar F0() const { return F0_; }
        
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
        )
        {
            this->calcThermFun(T, T0);
            return F_;
        }

        //- Return the F field interpolated to the element faces
        // Use an interpolation method instead of actually applying the
        // thermFunScheme over the face, to improve stability
        // (otherwise it may lead to mismatches betwen face and nodal values)
        virtual tmp<surfaceScalarField> Fface() const
        {
            return fvc::interpolate(F_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermFunSchemes 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
