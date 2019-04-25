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

#include "constitutiveEq.H"
#include <Eigen/Dense> // For eigen decomposition
#include "jacobi.H"    // Only required for jacobi decomposition

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constitutiveEq, 0);
    defineRunTimeSelectionTable(constitutiveEq, dictionary);
  
    template<>
    const char* NamedEnum
    <
        constitutiveEq::stabOptions,
        3
    >::names[] =
    {
        "none",
        "BSD",
        "coupling"
    };
  
    const NamedEnum<constitutiveEq::stabOptions, 3> constitutiveEq::stabOptionNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constitutiveEq::constitutiveEq
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& T
)
:
    name_(name),
    U_(U),
    phi_(phi),
    T_(T),
    rho_
    (
        IOobject
        (
            "rho",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zeroRho", dimDensity, 0)
    ),
    eta_
    (
        IOobject
        (
            "eta",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zeroEta", dimPressure*dimTime, 0)
    ),
    etaS_
    (
        IOobject
        (
            "etaS",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zeroEta", dimPressure*dimTime, 0)
    ),
    etaP_
    (
        IOobject
        (
            "etaP",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zeroEta", dimPressure*dimTime, 0)
    )
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> constitutiveEq::divTau
(   
    const volVectorField& U  
)
{
    volTensorField L = fvc::grad(U);
    if (isGNF())
    {   
        return  // accounts for variable eta
        (
            fvm::laplacian(eta()/rho(), U, "laplacian(eta,U)")
          + fvc::div(eta()/rho()*dev2(Foam::T(L)), "div(eta*dev2(T(gradU)))")
        ); 
    }
    else
    { 
        switch (stabOption_) 
        {
            case soNone :  // none - accounts for variable etaS or etaP
            return
            (
                fvc::div(tau()/rho(), "div(tau)")
              + fvm::laplacian(etaS()/rho(), U, "laplacian(eta,U)")
              + fvc::div(etaS()/rho()*dev2(Foam::T(L)), "div(eta*dev2(T(gradU)))")
            );

            case soBSD :  // BSD - accounts for variable etaS or etaP
            return
            (
                fvc::div(tau()/rho(), "div(tau)")
              + fvm::laplacian((etaP() + etaS())/rho(), U, "laplacian(eta,U)")
              + fvc::div(etaS()/rho()*dev2(Foam::T(L)), "div(eta*dev2(T(gradU)))")
              - fvc::laplacian(etaP()/rho(), U, "laplacian(etaP,U)")
            );

            case soCoupling :  // coupling - accounts for variable etaS or etaP
            return
            (
                fvc::div(tau()/rho(), "div(tau)")
              + fvm::laplacian((etaP() + etaS())/rho(), U, "laplacian(eta,U)")
              + fvc::div(etaS()/rho()*dev2(Foam::T(L)), "div(eta*dev2(T(gradU)))")
              - fvc::div(etaP()/rho()*L, "div(etaP*grad(U))")
              //- etaP() / rho() * fvc::div(L)
              //- (fvc::grad(etaP()/rho()) & L)
            );
         
            default:  // This will never happen
            return fvVectorMatrix(U, U.dimensions()); 
        }      
    }    
}

tmp<fvVectorMatrix> constitutiveEq::divTauS
(   
    const volVectorField& U,  
    const volScalarField& alpha
)
{      
    volTensorField L = fvc::grad(U);
    if (isGNF())  // accounts for variable eta
    {       
        return
        (
            fvm::laplacian(eta()*alpha, U, "laplacian(eta,U)")
          + fvc::div(eta()*alpha*dev2(Foam::T(L)), "div(eta*alpha*dev2(T(gradU)))")
        ); 
    }
    else
    {
        volTensorField L = fvc::grad(U);
        switch (stabOption_) 
        {
            case soNone :  // none - accounts for variable etaS or etaP
            return
            (        
                fvm::laplacian(etaS()*alpha, U, "laplacian(eta,U)")
              + fvc::div(etaS()*alpha*dev2(Foam::T(L)), "div(eta*alpha*dev2(T(gradU)))")
            );

            case soBSD :  // BSD - accounts for variable etaS or etaP
            return
            (
                fvm::laplacian((etaP() + etaS())*alpha, U, "laplacian(eta,U)")
              + fvc::div(etaS()*alpha*dev2(Foam::T(L)), "div(eta*alpha*dev2(T(gradU)))")
              - fvc::laplacian(etaP()*alpha, U, "laplacian(etaP,U)")
            );
          
            case soCoupling :  // coupling  - accounts for variable etaS or etaP
            return
            (        
                fvm::laplacian((etaP() + etaS())*alpha, U, "laplacian(eta,U)")
              + fvc::div(etaS()*alpha*dev2(Foam::T(L)), "div(eta*alpha*dev2(T(gradU)))")
              - fvc::div(etaP()*alpha*L, "div(etaP*grad(U))")
            ); 
         
            default:  // This will never happen
            return fvVectorMatrix(U, U.dimensions());  
        }      
    }      
}

void constitutiveEq::decomposeGradU
(
    const volTensorField& M,
    const volTensorField& eigVals, 
    const volTensorField& eigVecs,
    volTensorField& omega, 
    volTensorField& B
)
{
    forAll(M, cellI)
    {
        const tensor& eigValsR = eigVals[cellI];
        const tensor& MR = M[cellI];
        tensor& omegaR = omega[cellI];
            
        B[cellI].xx()=MR.xx();
        B[cellI].yy()=MR.yy(); 
        B[cellI].zz()=MR.zz(); 
             
        omegaR.xy() = ( eigValsR.yy()*MR.xy() + eigValsR.xx()*MR.yx() )
            / ( eigValsR.yy() - eigValsR.xx() + 1e-16 );
        omegaR.xz() = ( eigValsR.zz()*MR.xz() + eigValsR.xx()*MR.zx() )
            / ( eigValsR.zz() - eigValsR.xx() + 1e-16 );
        omegaR.yz() = ( eigValsR.zz()*MR.yz() + eigValsR.yy()*MR.zy() )
            / ( eigValsR.zz() - eigValsR.yy() + 1e-16 );

        omegaR.yx() = -omegaR.xy();
        omegaR.zx() = -omegaR.xz();
        omegaR.zy() = -omegaR.yz(); 
    }

    omega = ( eigVecs & omega & eigVecs.T() );

    B = ( eigVecs & B & eigVecs.T() );
}

void constitutiveEq::calcEig
(
    const volSymmTensorField& theta,
    volTensorField& vals,
    volTensorField& vecs
)
{
    // Eigen decomposition using a QR algorithm of Eigen library 
    Eigen::Matrix3d theta_eig(Eigen::Matrix3d::Zero(3,3));
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;

    forAll(theta, cellI)
    {
        // Transfer theta from OF to Eigen        
        const symmTensor& thetaR = theta[cellI];  
           
        theta_eig(0,0)=thetaR.xx();
        theta_eig(1,1)=thetaR.yy();
        theta_eig(2,2)=thetaR.zz();

        theta_eig(0,1)=thetaR.xy();
        theta_eig(1,0)=thetaR.xy();

        theta_eig(0,2)=thetaR.xz();
        theta_eig(2,0)=thetaR.xz();

        theta_eig(1,2)=thetaR.yz();
        theta_eig(2,1)=thetaR.yz();

        // Compute eigenvalues/vectors in Eigen
        eigSol.compute(theta_eig);
        Eigen::Vector3d eival = eigSol.eigenvalues();
        Eigen::Matrix3d eivect = eigSol.eigenvectors();

        // Transfer eigenvalues/vectors from Eigen to OF 
        tensor& vecsR = vecs[cellI];  

        vecsR.xx()=eivect(0,0);
        vecsR.yx()=eivect(1,0);
        vecsR.zx()=eivect(2,0);

        vecsR.xy()=eivect(0,1);
        vecsR.yy()=eivect(1,1);      
        vecsR.zy()=eivect(2,1);      

        vecsR.xz()=eivect(0,2);  
        vecsR.yz()=eivect(1,2);
        vecsR.zz()=eivect(2,2);

        vals[cellI] *= 0.;
        vals[cellI].xx()=Foam::exp(eival(0));
        vals[cellI].yy()=Foam::exp(eival(1));
        vals[cellI].zz()=Foam::exp(eival(2));
    }
/*
 // Eigen decomposition using the iterative jacobi algorithm 

    forAll(theta, cellI)
    {
        int N=3;
        int NROT=0;
        jacobi(theta[cellI], N, vals[cellI], vecs[cellI], NROT);
    }
*/
}

void constitutiveEq::checkForStab(const dictionary& dict)
{
    stabOption_ = stabOptionNames_.read(dict.lookup("stabilization"));
}

tmp<volSymmTensorField> constitutiveEq::tauTotal()
{
    volTensorField L = fvc::grad(U_);
     
    if (isGNF())
    { 
        //return eta()*symm(L+L.T());  //original, incompressible
        //return eta()*2*symm(L);  //incompressible
        return eta()*dev(2*symm(L));  //compressible
    }
    else
    {
        //return tau() + etaS()*symm(L+L.T());  //original, incompressible
        //return tau() + etaS()*2*symm(L);  //incompressible
        return tau() + etaS()*dev(2*symm(L));  //compressible
    }
}

tmp<volSymmTensorField> constitutiveEq::tauN()
{
    volTensorField L = fvc::grad(U_);
     
    if (isGNF())
    { 
        //return eta()*2*symm(L);  //incompressible
        return eta()*dev(2*symm(L));  //compressible
    }
    else
    {
        //return etaS()*2*symm(L);  //incompressible
        return etaS()*dev(2*symm(L));  //compressible
    }
}

} //End namespace

// ************************************************************************* //
