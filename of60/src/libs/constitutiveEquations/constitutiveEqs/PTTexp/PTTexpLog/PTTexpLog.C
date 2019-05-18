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

#include "PTTexpLog.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(PTTexpLog, 0);
    addToRunTimeSelectionTable(constitutiveEq, PTTexpLog, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::PTTexpLog::PTTexpLog
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
    theta_
    (
        IOobject
        (
            "theta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    eigVals_
    (
        IOobject
        (
            "eigVals" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),
    eigVecs_
    (
        IOobject
        (
            "eigVecs" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
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
        Info<< "Neglecting fluid structure dependence on temperature" 
            << " in the Constitutive Eq" << endl;
        calcDotTHfun_ = 0;
    }
    Info<< "calcDotTHfun_ = " << calcDotTHfun_ << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::PTTexpLog::dotTHfun()
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
            elastEnergDiss_/T()*DTDt - fvc::div(phi())
        )
    );
}

void Foam::constitutiveEqs::PTTexpLog::correct()
{
    dimensionedSymmTensor Itensor
    (
        "Identity", dimensionSet(0, 0, 0, 0, 0, 0, 0), symmTensor::I
    );

    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    dimensionedScalar c1( "zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.);
    volTensorField B = c1 * eigVecs_;
    volTensorField omega = B;
    volTensorField M = (eigVecs_.T() & ( L.T() - zeta_*symm(L) ) & eigVecs_);
    decomposeGradU (M, eigVals_, eigVecs_, omega, B);

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

    // Solve the constitutive Eq in theta = log(c)
    fvSymmTensorMatrix thetaEqn
    (
        fvm::ddt(theta_)
      + fvm::div(phi(), theta_)
      ==
        symm
        (
            (omega & theta_) - (theta_ & omega) + 2.0 * B
          + Foam::exp(epsilon_ / (1.0 - zeta_)
          * (tr((eigVecs_ & eigVals_ & eigVecs_.T())) - 3.0)) / lambda_
          * ( eigVecs_ & (inv(eigVals_) - Itensor) & eigVecs_.T() )
          + calcDotTHfun_ * dotTHfun() * Itensor  // Thermal dependency dotT*H
          // The following term appears when the constitutive equation is formulated
          // in terms of deviatoric tau tensor instead of the Cauchy sigma tensor,
          // when expanding the upper-convected derivative of sigma.
          // It yields zero if the thermal functions of lambda and etaP are the same.
          + (dotLambda / lambda_ - dotEtaP / etaP()) * Itensor
        )
    );

    thetaEqn.relax();
    thetaEqn.solve();
    // Diagonalization of theta
    calcEig(theta_, eigVals_, eigVecs_);
    // Convert from theta to tau
    tau_ = (etaP() / lambda_ / (1.0 - zeta_))
        * symm((eigVecs_ & eigVals_ & eigVecs_.T()) - Itensor);
    tau_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::constitutiveEqs::PTTexpLog::energyExtraTerms()
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
                * Foam::exp(epsilon_*lambda_/etaP()*tr(tau_))
                * tr(tau_)/(2.0*lambda_)
            )
        )
    );
}

// ************************************************************************* //
