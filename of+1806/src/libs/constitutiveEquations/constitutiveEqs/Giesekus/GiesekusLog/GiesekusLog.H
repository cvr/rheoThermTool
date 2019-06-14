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
    GiesekusLog

Description
    GiesekusLog viscoelastic fluid model with inclusion of thermal effects.
    
    Based from:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves.
    
    This class is part of rheoThermTool.

SourceFiles
    GiesekusLog.C

\*---------------------------------------------------------------------------*/

#ifndef GiesekusLog_H
#define GiesekusLog_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class GiesekusLog Declaration
\*---------------------------------------------------------------------------*/

class GiesekusLog
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Logarithm of the conformation tensor
        volSymmTensorField theta_;

        //- Eigenvalues of theta diagonalization (values on the tensor diagonal)
        volTensorField eigVals_;

        //- Eigenvectors of theta diagonalization
        volTensorField eigVecs_;

        // Model constants

            //- elasticEnergyDissipation factor
            dimensionedScalar elastEnergDiss_;

            //- Mobility factor
            dimensionedScalar alpha_;

            //- Relaxation time
            volScalarField lambda_;

         // Model switches

            bool dotLambdaSwitch_;
            
            bool dotEtaPSwitch_;

            bool dotTHfunSwitch_;
            
            dimensionedScalar calcDotTHfun_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GiesekusLog(const GiesekusLog&);

        //- Disallow default bitwise assignment
        void operator=(const GiesekusLog&);

public:

    //- Runtime type information
    TypeName("GiesekusLog");

    // Constructors

        //- Construct from components
        GiesekusLog
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~GiesekusLog() {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const { return tau_; }

        //- Return the elastic relaxation time
        //  Return reference for manipulation outside constitutiveEq Class
        virtual tmp<volScalarField> lambda() const { return lambda_; }
        //virtual volScalarField* lambdaPtr() { return &lambda_; }
        virtual volScalarField& lambdaRef() { return lambda_; }

        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const { return false; }

        //- Thermal dependency function for the constitutive equation
        virtual tmp<volScalarField> dotTHfun();

        //- Correct the viscoelastic stress
        virtual void correct();

        //- Extra viscoelastic dissipation terms for the energy equation
        virtual tmp<volScalarField> energyExtraTerms();


    // Null ptr and refs for manipulation outside constitutiveEq Class
        
        //- GNF
        //virtual volScalarField& etaRef() { return this->eta(); }
        virtual volScalarField* etaPtr() { return nullptr; }
        virtual volScalarField& etaRef() { return *etaPtr(); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
