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
    constitutiveModel

Description
    Base class which can be directly called by applications to use any of 
    the models derived from class constitutivEq.

    Based from:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves. 
    
    This class is part of rheoThermTool.

SourceFiles
    constitutiveModel.C

\*---------------------------------------------------------------------------*/

#ifndef constitutiveModel_H
#define constitutiveModel_H

#include "constitutiveEq.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constitutiveModel Declaration
\*---------------------------------------------------------------------------*/

class constitutiveModel
:
    public IOdictionary
{
private:

    // Private data

        //- Constitutive Eq
        autoPtr<constitutiveEq> eqPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveModel(const constitutiveModel&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveModel&);

public:

    //- Runtime type information
    TypeName("constitutiveModel");


    // Constructors

        //- Construct from components
        constitutiveModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T
        );

    // Destructor

        virtual ~constitutiveModel()
        {}


    // Member Functions

        //- Return the constitutiveEq model name
        const word& eqTypeName() const { return eqPtr_->type(); }

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() { return eqPtr_->tau(); }
        
        //- Return the total stress tensor (solvent + polymeric)
        virtual tmp<volSymmTensorField> tauTotal() { return eqPtr_->tauTotal(); }
        
        //- Return the solvent stress tensor (newtonian)
        virtual tmp<volSymmTensorField> tauN() { return eqPtr_->tauN(); }
        
        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U)
        {
            return eqPtr_->divTau(U);
        }
        
        //- Return the density
        virtual tmp<volScalarField> rho() const { return eqPtr_->rho(); }
        //virtual volScalarField* rhoPtr() { return eqPtr_->rhoPtr(); }
        virtual volScalarField& rhoRef() { return eqPtr_->rhoRef(); }
        
        //- Return the viscosity field for GNF models
        virtual tmp<volScalarField> eta() const { return eqPtr_->eta(); }
        //virtual volScalarField* etaPtr() { return eqPtr_->etaPtr(); }
        virtual volScalarField& etaRef() { return eqPtr_->etaRef(); }

        //- Return the solvent viscosity for VE models
        virtual tmp<volScalarField> etaS() const { return eqPtr_->etaS(); }
        //virtual volScalarField* etaSPtr() { return eqPtr_->etaSPtr(); }
        virtual volScalarField& etaSRef() { return eqPtr_->etaSRef(); }

        //- Return the polymeric viscosity for VE models
        virtual tmp<volScalarField> etaP() const { return eqPtr_->etaP(); }
        //virtual volScalarField* etaPPtr() { return eqPtr_->etaPPtr(); }
        virtual volScalarField& etaPRef() { return eqPtr_->etaPRef(); }

        //- Is GNF or VE
        virtual bool isGNF() const { return eqPtr_->isGNF(); }
       
        //- Correct the viscoelastic stress
        virtual void correct() { eqPtr_->correct(); }

	//- Extra viscoelastic dissipation terms for the energy equation
	virtual tmp<volScalarField> energyExtraTerms()
	{
	    return eqPtr_->energyExtraTerms();
	}

        //- Read constitutiveProperties dictionary
        virtual bool read();


    // Null ptr and refs for manipulation outside constitutiveEq Class
        
        //- GNF fields
        virtual volScalarField* kPtr() { return eqPtr_->kPtr(); }
        virtual volScalarField& kRef() { return eqPtr_->kRef(); }
        virtual volScalarField* etaMinPtr() { return eqPtr_->etaMinPtr(); }
        virtual volScalarField& etaMinRef() { return eqPtr_->etaMinRef(); }
        virtual volScalarField* etaMaxPtr() { return eqPtr_->etaMaxPtr(); }
        virtual volScalarField& etaMaxRef() { return eqPtr_->etaMaxRef(); }

        //- VE fields
        virtual volScalarField* lambdaPtr() { return eqPtr_->lambdaPtr(); }
        virtual volScalarField& lambdaRef() { return eqPtr_->lambdaRef(); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
