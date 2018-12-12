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
    CarreauYasuda

Description
    Carreau-Yasuda flow model, based on:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves.

SourceFiles
    CarreauYasuda.C

\*---------------------------------------------------------------------------*/

#ifndef CarreauYasuda_H
#define CarreauYasuda_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class CarreauYasuda Declaration
\*---------------------------------------------------------------------------*/

class CarreauYasuda
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        // Model constants

            //- Dimensionaless parameters
            dimensionedScalar a_;
            dimensionedScalar n_;

            //- Characteristic time-scale (3d if function of T)
            volScalarField k_;

            //- Zero strain rate viscosity (eta0, 3d if function of T) 
            volScalarField etaMin_;

            //- Infinite strain rate viscosity (etaInf, 3d if function of T)
            volScalarField etaMax_;


    // Private Member Functions
        //- Disallow default bitwise copy construct
        CarreauYasuda(const CarreauYasuda&);

        //- Disallow default bitwise assignment
        void operator=(const CarreauYasuda&);

public:

    //- Runtime type information
    TypeName("CarreauYasuda");

    // Constructors

        //- Construct from components
        CarreauYasuda
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~CarreauYasuda()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const { return tau_; }

        //- Return consistency index
        //  Return reference for manipulation outside constitutiveEq Class
        virtual tmp<volScalarField> k() const { return k_; }
        //virtual volScalarField* kPtr() { return &k; }
        virtual volScalarField& kRef() { return k_; }

        //- Return Lower bound for viscosity
        //  Return reference for manipulation outside constitutiveEq Class
        virtual tmp<volScalarField> etaMin() const { return etaMin_; }
        //virtual volScalarField* etaMinPtr() { return &etaMin_; }
        virtual volScalarField& etaMinRef() { return etaMin_; }

        //- Upper bound for viscosity
        //  Return reference for manipulation outside constitutiveEq Class
        virtual tmp<volScalarField> etaMax() const { return etaMax_; }
        //virtual volScalarField* etaMaxPtr() { return &etaMax_; }
        virtual volScalarField& etaMaxRef() { return etaMax_; }

        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const { return true; }

        //- Correct the viscoelastic stress
        virtual void correct();

        //- Extra viscoelastic dissipation terms for the energy equation
        virtual tmp<volScalarField> energyExtraTerms();


    // Null ptr and refs for manipulation outside constitutiveEq Class
        
        //- VE fields
        //virtual volScalarField& etaSRef() { return this->etaS(); }
        virtual volScalarField* etaSPtr() { return nullptr; }
        virtual volScalarField& etaSRef() { return *etaSPtr(); }
        //virtual volScalarField& etaPRef() { return this->etaP(); }
        virtual volScalarField* etaPPtr() { return nullptr; }
        virtual volScalarField& etaPRef() { return *etaPPtr(); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
