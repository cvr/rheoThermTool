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
    multiMode

Description
    Software for combining several viscoelastic fluid models together,
    as a single fluid composed of several sets of solvents and polymer chains.
    It has support for the thermal effects included in each model,
    but may lack full theoretical support.
    Original code based from:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves.

SourceFiles
    multiMode.C

\*---------------------------------------------------------------------------*/

#ifndef multiMode_H
#define multiMode_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class multiMode Declaration
\*---------------------------------------------------------------------------*/

class multiMode
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        mutable volSymmTensorField tau_;

        //- List of models
        PtrList<constitutiveEq> models_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMode(const multiMode&);

        //- Disallow default bitwise assignment
        void operator=(const multiMode&);

public:

    //- Runtime type information
    TypeName("multiMode");

    // Constructors

        //- Construct from components
        multiMode
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMode()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const;
        
        //- Return the density as the average density field
        //  considering each model composing the multimode.
        virtual void rhoMulti();

        //- Return the etaS as the sum of the etaS field
        //  considering each model composing the multimode.
        //  Note v3.0: etaS is needed for tauTotal()
        virtual void etaSMulti();
    
        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U);
        
        //- Return the solvent extra-stress contribution + stabilizing terms for phase alpha
        virtual tmp<fvVectorMatrix> divTauS(const volVectorField& U, const volScalarField& alpha);
        
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const { return false; }

        //- Correct the viscoelastic stress
        virtual void correct();

        //- Extra viscoelastic dissipation terms for the energy equation
        virtual tmp<volScalarField> energyExtraTerms() ;


    // Null ptr and refs for manipulation outside constitutiveEq Class
        
        //- GNF fields
        //virtual volScalarField& etaRef() { return this->eta(); }
        virtual volScalarField* etaPtr() { return nullptr; }
        virtual volScalarField& etaRef() { return *etaPtr(); }

        //- VE fields
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
