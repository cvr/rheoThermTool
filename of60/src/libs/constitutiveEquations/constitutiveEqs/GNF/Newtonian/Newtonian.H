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
    Newtonian

Description
    Newtonian flow model.
    
    Based from:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves.
    
    This class is part of rheoThermTool.

SourceFiles
    Newtonian.C

\*---------------------------------------------------------------------------*/

#ifndef Newtonian_H
#define Newtonian_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class Newtonian
:
    public constitutiveEq
{
    // Private data
        
        //- Transported viscoelastic stress
        volSymmTensorField tau_;


    // Private Member Functions
        //- Disallow default bitwise copy construct
        Newtonian(const Newtonian&);

        //- Disallow default bitwise assignment
        void operator=(const Newtonian&);

public:

    //- Runtime type information
    TypeName("Newtonian");

    // Constructors

        //- Construct from components
        Newtonian
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~Newtonian()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const { return tau_; }
        
        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U) const;
        
        //- Return the solvent extra-stress contribution + stabilizing terms
        //  for phase alpha
        virtual tmp<fvVectorMatrix> divTauS
        (
            const volVectorField& U,
            const volScalarField& alpha
        ) const;

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
