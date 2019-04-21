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

\*---------------------------------------------------------------------------*/

#include "PTTlinear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(PTTlinear, 0);
    addToRunTimeSelectionTable(constitutiveEq, PTTlinear, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::PTTlinear::PTTlinear
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& T,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi, T),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    elastEnergDiss_(dict.lookup("elastEnergDiss")),
    epsilon_(dict.lookup("epsilon")),
    zeta_(dict.lookup("zeta")),
    lambda_
    (
        IOobject
        (
            "lambda" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("lambda"))
    ),
    dotLambdaSwitch_(dict.lookupOrDefault<Switch>("dotLambdaSwitch", true)),
    dotEtaPSwitch_(dict.lookupOrDefault<Switch>("dotEtaPSwitch", true)),
    dotTHfunSwitch_(dict.lookupOrDefault<Switch>("dotTHfunSwitch", true)),
    calcDotTHfun_("calcDotTHfun", dimless, 1)
{
    checkForStab(dict);

    rhoRef().rename("rho" + name);
    rhoRef().readOpt() = IOobject::READ_IF_PRESENT;
    rhoRef().writeOpt() = IOobject::AUTO_WRITE;
    rhoRef() = volScalarField
    (
        IOobject
        (
            "rho" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("rho"))
    );

    etaSRef().rename("etaS" + name);
    etaSRef().readOpt() = IOobject::READ_IF_PRESENT;
    etaSRef().writeOpt() = IOobject::AUTO_WRITE;
    etaSRef() = volScalarField
    (
        IOobject
        (
            "etaS" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("etaS"))
    );

    etaPRef().rename("etaP" + name);
    etaPRef().readOpt() = IOobject::READ_IF_PRESENT;
    etaPRef().writeOpt() = IOobject::AUTO_WRITE;
    etaPRef() = volScalarField
    (
        IOobject
        (
            "etaP" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dict.lookup("etaP"))
    );

    if (!dotLambdaSwitch_)
    {
        Info<< "Neglecting influence of lambda total derivative" << endl;
    }
    
    if (!dotEtaPSwitch_)
    {
        Info<< "Neglecting influence of etaP total derivative" << endl;
    }

    if (!dotTHfunSwitch_)
    {
        Info<< "Neglecting fluid structure dependence on temperature" << endl;
        calcDotTHfun_ = 0;
    }
    Info<< "calcDotTHfun_ = " << calcDotTHfun_ << endl;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::PTTlinear::dotTHfun()
{
    //- Get dimensionless ln(T) field to compute DlnT/Dt
    //  WARNING: This may yield problems if the ln(T) field is not
    //  a properly build with an IOobject due to the time derivative!
    //dimensionedScalar Tone (word::null, T().dimensions(), 1); 
    //volScalarField lnT ("lnT", Foam::log(T/Tone));

    // Compute DT/Dt
    //volScalarField DTDt = fvc::ddt(T())
    //    + fvc::div(phi(),T()) - T()*fvc::div(phi());
    volScalarField DTDt = fvc::DDt(phi(), T());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "dotTHfun",
            calcDotTHfun_ * ( elastEnergDiss_/T()*DTDt - fvc::div(phi()) )
        )
    );
}

void Foam::constitutiveEqs::PTTlinear::correct()
{
    dimensionedSymmTensor Itensor
    ( 
        "Identity", dimensionSet(0, 0, 0, 0, 0, 0, 0), symmTensor::I
    );

    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    // Compute dotLambda = Dlambda/Dt
    //volScalarField dotLambda = fvc::ddt(lambda_)
    //    + fvc::div(phi(),lambda_) - lambda_*fvc::div(phi());
    volScalarField dotLambda = fvc::DDt(phi(), lambda_);
    if (!dotLambdaSwitch_) { dotLambda *= 0; }
        
    // Compute dotEtaP = DetaP/Dt
    //volScalarField dotEtaP = fvc::ddt(etaPRef())
    //    + fvc::div(phi(),etaPRef()) - etaPRef()*fvc::div(phi());
    volScalarField dotEtaP = fvc::DDt(phi(), etaPRef());
    if (!dotEtaPSwitch_) { dotEtaP *= 0; }

    //Info<< "DEBUG max(|dotLambda|) = " << max(mag(dotLambda)) << endl;
    //Info<< "DEBUG max(|dotEtaP|) = " << max(mag(dotEtaP)) << endl;
    //Info<< "DEBUG max(|dotTHfun|) = " << max(mag(dotTHfun())) << endl;

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        twoSymm(C)
      - fvm::Sp(epsilon_/etaP()*tr(tau_) + 1.0/lambda_, tau_)
      //- 0.5*zeta_*(symm(tau_ & twoD) + symm(twoD & tau_))  // why? ...
      - zeta_*(symm(tau_ & twoD))  // this must yield the same as above!
      //+ etaP()/lambda_*twoD  // incorrect, check below
      + (1.0 - zeta_)*etaP()/lambda_*twoD  // slip factor should also be here
      // Thermal dependency dotT*Ht
      - fvm::SuSp(-dotTHfun(), tau_)  // Thermal dependency dotT*H*tau_
      + dotTHfun()*etaP()/lambda_*Itensor  // Thermal dependency dotT*H*G*I
      // The following term appears when the constitutive equation is formulated
      // in terms of deviatoric tau tensor instead of the Cauchy sigma tensor,
      // when expanding the upper-convected derivative of sigma.
      // It yields zero if the thermal functions of lambda and etaP are the same.
      + (etaP()/lambda_*dotLambda - dotEtaP)/lambda_*Itensor
    );

    tauEqn.relax();
    tauEqn.solve();
}

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::PTTlinear::energyExtraTerms()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "energyExtraTerms",
            (
                (zeta_ + elastEnergDiss_ - zeta_*elastEnergDiss_)
                * (
                    //(tau_ && symm(fvc::grad(U())))  // can be simplified...
                    (tau_ && fvc::grad(U()))  // T:S==T:grad(U) for symmetric T
                  + etaP()/lambda_*fvc::div(phi())
                  )
              + (1.0 - elastEnergDiss_)
                * (epsilon_*lambda_/etaP()*tr(tau_) + 1.0)
                * tr(tau_)/(2.0*lambda_)
            )
        )
    );
}

// ************************************************************************* //
