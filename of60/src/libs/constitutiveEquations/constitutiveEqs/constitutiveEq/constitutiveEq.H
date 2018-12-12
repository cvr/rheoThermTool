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
    constitutiveEq

Description
    Abstract base class from which each specific GNF or viscoelastic model
    is derived.

    Based on:
    - libviscoelasticTransportModels library written by Jovani L. Favero.
    - rheoTool library written by Francisco Pimenta and Manuel A. Alves. 

SourceFiles
    constitutiveEq.C
    constitutiveEqNew.C

\*---------------------------------------------------------------------------*/

#ifndef constitutiveEq_H
#define constitutiveEq_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "NamedEnum.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constitutiveEq Declaration
\*---------------------------------------------------------------------------*/

class constitutiveEq
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;
        
        //- Reference to temperature field
        const volScalarField& T_;

        //- Density (3D if function of EOS)
        volScalarField rho_;

        //- GNF viscosity field (3D if function of T)
        //  (default behavior is a null field - for VE models)  
        volScalarField eta_;

        //- VE solvent viscosity (3D if function of T)
        //  (default behavior is a null field - for GNF models)  
        volScalarField etaS_;

        //- VE zero shear rate polymer viscosity (3D if function of T)
        //  (default behavior is a null field - for GNF models)  
        volScalarField etaP_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveEq(const constitutiveEq&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveEq&);
        
public:

    //- Stabilization options
    enum stabOptions
    {
      soNone,
      soBSD,
      soCoupling
    };
        
protected:

    // Protected Members
   
        //- Return B and Sigma matrices after decomposition of grad(U).T()
        //  for the LCM
        void decomposeGradU
        (
            const volTensorField& M,
            const volTensorField& eigVals, 
            const volTensorField& eigVecs,
            volTensorField& sigma, 
            volTensorField& B
        );

        //- Return eigenvectors/values of theta for the LCM
        void calcEig
        (
            const volSymmTensorField& theta,
            volTensorField& vals,
            volTensorField& vecs
        );
        
        //- Return the strain rate magnitude for GNF models
        inline tmp<volScalarField> strainRate()
        {
            return ( sqrt(2.0)*mag(symm(fvc::grad(U_))) );
        }
        
        //- Check which stabilization method to use
        void checkForStab(const dictionary& dict);
      
        static const NamedEnum<stabOptions, 3> stabOptionNames_;
        stabOptions stabOption_;

public:

    //- Runtime type information
    TypeName("constitutiveEq");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            constitutiveEq,
            dictionary,
            (
                const word& name,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const volScalarField& T,
                const dictionary& dict
            ),
            (name, U, phi, T, dict)
        );


    // Constructors

        //- Construct from components
        constitutiveEq
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<constitutiveEq> New
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~constitutiveEq()
        {}

    // Member Functions

        //- Return name
        const word& name() const { return name_; }

        //- Return the velocity field
        const volVectorField& U() const { return U_; }

        //- Return the face flux field
        const surfaceScalarField& phi() const { return phi_; }

        //- Return the temperature field
        const volScalarField& T() const { return T_; }

    // The following functions are intended to be used inside divTau() and
    // divTauS() functions, which should be general for almost all models
    // (exception: multiMode). Any other use should be carefully analyzed,
    // with special attention to the mentioned exceptions.
    //
    // The references are also returned (*Ref functions) to allow
    // for the field to be edited outside the constitutiveEq class

        //- Return the density
        virtual tmp<volScalarField> rho() const { return rho_; }
        //virtual volScalarField* rhoPtr() { return &rho_; }
        virtual volScalarField& rhoRef() { return rho_; }

        //- Return the viscosity field for GNF models
        //  (default behavior is a null field - for VE models)  
	virtual tmp<volScalarField> eta() const { return eta_; }
	//virtual volScalarField* etaPtr() { return &eta_; }
	virtual volScalarField& etaRef() { return eta_; }

        //- Return the solvent viscosity for VE models
        virtual tmp<volScalarField> etaS() const { return etaS_; }
	//virtual volScalarField* etaSPtr() { return &etaS_; }
	virtual volScalarField& etaSRef() { return etaS_; }

        //- Return the polymeric viscosity for VE models
        //  Does not consider possible variation of etaP with shear rate
        //  (there is no need for it where is used)
	virtual tmp<volScalarField> etaP() const { return etaP_; }
	//virtual volScalarField* etaPPtr() { return &etaP_; }
        virtual volScalarField& etaPRef() { return etaP_; }


    // Other functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;

        //- Return the coupling term for the momentum equation (since it is
        //  virtual non-pure, it does not need to be redefined in derived
        //  classes, although it can if any model needs a
        //  different implementation)
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U);
        
        //- Return the solvent extra-stress contribution + stabilizing terms
        //  for phase alpha (since it is virtual non-pure, it does not need to
        //  be redefined in derived classes, although it can if any model needs
        //  a different implementation)
        virtual tmp<fvVectorMatrix> divTauS
        (
	    const volVectorField& U,
	    const volScalarField& alpha
	);
        
        //- Return the total stress tensor (solvent + polymeric)
        //  using alpha1 as weight
        virtual tmp<volSymmTensorField> tauTotal();
        
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const = 0;
        
        //- Thermal dependency function for the constitutive equation
        //virtual tmp<volScalarField> dotTHfun() = 0;

        //- Correct the viscoelastic stress
        virtual void correct() = 0;

        //- Extra viscoelastic dissipation terms for the energy equation
        virtual tmp<volScalarField> energyExtraTerms() = 0;


    // Null ptr and refs for manipulation outside constitutiveEq Class
        
        //- GNF fields
        virtual volScalarField* kPtr() { return nullptr; }
        virtual volScalarField& kRef() { return *kPtr(); }
        virtual volScalarField* etaMinPtr() { return nullptr; }
        virtual volScalarField& etaMinRef() { return *etaMinPtr(); }
        virtual volScalarField* etaMaxPtr() { return nullptr; }
        virtual volScalarField& etaMaxRef() { return *etaMaxPtr(); }
        
        //- VE fields
        virtual volScalarField* lambdaPtr() { return nullptr; }
        virtual volScalarField& lambdaRef() { return *lambdaPtr(); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
