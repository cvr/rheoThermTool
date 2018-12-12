/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::constitutiveTwoPhaseMixture

Description
    A two-phase GNF or VE transportModel

SourceFiles
    constitutiveTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef constitutiveTwoPhaseMixture_H
#define constitutiveTwoPhaseMixture_H

#include "constitutiveEq.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class constitutiveTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class constitutiveTwoPhaseMixture
:
    public IOdictionary,
    public twoPhaseMixture
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveTwoPhaseMixture(const constitutiveTwoPhaseMixture&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveTwoPhaseMixture&);


protected:

    // Protected data

        autoPtr<constitutiveEq> phase1_;
        autoPtr<constitutiveEq> phase2_;
        
        /*
        volScalarField rho1_;
        volScalarField rho2_;
        */

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        const volScalarField& T_;
        
        // Polymeric extra-stress
        volSymmTensorField tauMF_;


public:

    TypeName("constitutiveTwoPhaseMixture");


    // Constructors

        //- Construct from components
        constitutiveTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& T
        );


    //- Destructor
    virtual ~constitutiveTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 constitutiveModel
        const constitutiveEq& phase1() const { return phase1_(); }

        //- Return const-access to phase2 constitutiveModel
        const constitutiveEq& phase2() const { return phase2_(); }
        
        /*
        //- Return const-access to phase1 density
        const tmp<volScalarField> rho1() const { return rho1_; }
        
        //- Return const-access to phase2 density
        const tmp<volScalarField> rho2() const { return rho2_; };
        */
        
        
        //- Return true if phase1 is GNF (non-elastic)
        virtual bool isGNF1() const { return phase1_->isGNF(); }

        //- Return true if phase2 is GNF (non-elastic)
        virtual bool isGNF2() const { return phase2_->isGNF(); }


        //- Return phase1 and phase2 density (access, pre, ref)
        virtual tmp<volScalarField> rho1() const { return phase1_->rho(); }
        //virtual volScalarField* rho1Ptr() { return phase1_->rhoPtr(); }
        virtual volScalarField& rho1Ref() { return phase1_->rhoRef(); }

        virtual tmp<volScalarField> rho2() const { return phase2_->rho(); }
        //virtual volScalarField* rho2Ptr() { return phase2_->rhoPtr(); }
        virtual volScalarField& rho2Ref() { return phase2_->rhoRef(); }


        //- Return phase1 and phase2 eta (access, pre, ref)
        virtual tmp<volScalarField> eta1() const { return phase1_->eta(); }
        //virtual volScalarField* eta1Ptr() { return phase1_->etaPtr(); }
        virtual volScalarField& eta1Ref() { return phase1_->etaRef(); }
        
        virtual tmp<volScalarField> eta2() const { return phase2_->eta(); }
        //virtual volScalarField* eta2Ptr() { return phase2_->etaPtr(); }
        virtual volScalarField& eta2Ref() { return phase2_->etaRef(); }


        //- Return phase1 and phase2 etaS (access, pre, ref)
        virtual tmp<volScalarField> etaS1() const { return phase1_->etaS(); }
        //virtual volScalarField* etaS1Ptr() { return phase1_->etaSPtr(); }
        virtual volScalarField& etaS1Ref() { return phase1_->etaSRef(); }
        
        virtual tmp<volScalarField> etaS2() const { return phase2_->etaS(); }
        //virtual volScalarField* etaS2Ptr() { return phase2_->etaSPtr(); }
        virtual volScalarField& etaS2Ref() { return phase2_->etaSRef(); }

        
        //- Return phase1 and phase2 etaP (access, pre, ref)
        virtual tmp<volScalarField> etaP1() const { return phase1_->etaP(); }
        //virtual volScalarField* etaP1Ptr() { return phase1_->etaPPtr(); }
        virtual volScalarField& etaP1Ref() { return phase1_->etaPRef(); }
        
        virtual tmp<volScalarField> etaP2() const { return phase2_->etaP(); }
        //virtual volScalarField* etaP2Ptr() { return phase2_->etaPPtr(); }
        virtual volScalarField& etaP2Ref() { return phase2_->etaPRef(); }
        

        //- Update the constitutiveEqs
        void correct()
        {
            phase1_->correct();
            phase2_->correct();
        }
        
        //- Polymeric extra-stress tensor
        const tmp<volSymmTensorField> tauMF() const { return tauMF_; }
      
        //- Return the multiphase coupling term for the momentum equation
        tmp<fvVectorMatrix> divTauMF(const volVectorField& U);
        
        //- Total extra-stress tensor (polymeric + solvent)
        tmp<volSymmTensorField> tauTotalMF();
        
        //- Read base transportProperties dictionary
        virtual bool read();


    // Null ptr and refs for manipulation outside constitutiveEq Class

        //- Return phase1 and phase2 lambda (access, pre, ref)
        virtual volScalarField* lambda1Ptr() { return phase1_->lambdaPtr(); }
        virtual volScalarField& lambda1Ref() { return phase1_->lambdaRef(); }
        
        virtual volScalarField* lambda2Ptr() { return phase2_->lambdaPtr(); }
        virtual volScalarField& lambda2Ref() { return phase2_->lambdaRef(); }
        
        
        //- Return phase1 and phase2 k (access, pre, ref)
        virtual volScalarField* k1Ptr() { return phase1_->kPtr(); }
        virtual volScalarField& k1Ref() { return phase1_->kRef(); }
        
        virtual volScalarField* k2Ptr() { return phase2_->kPtr(); }
        virtual volScalarField& k2Ref() { return phase2_->kRef(); }


        //- Return phase1 and phase2 etaMin (access, pre, ref)
        virtual volScalarField* etaMin1Ptr() { return phase1_->etaMinPtr(); }
        virtual volScalarField& etaMin1Ref() { return phase1_->etaMinRef(); }
        
        virtual volScalarField* etaMin2Ptr() { return phase2_->etaMinPtr(); }
        virtual volScalarField& etaMin2Ref() { return phase2_->etaMinRef(); }


        //- Return phase1 and phase2 etaMax (access, pre, ref)
        virtual volScalarField* etaMax1Ptr() { return phase1_->etaMaxPtr(); }
        virtual volScalarField& etaMax1Ref() { return phase1_->etaMaxRef(); }
        
        virtual volScalarField* etaMax2Ptr() { return phase2_->etaMaxPtr(); }
        virtual volScalarField& etaMax2Ref() { return phase2_->etaMaxRef(); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
