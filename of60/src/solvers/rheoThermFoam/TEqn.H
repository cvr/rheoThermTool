// Solve the energy equation

//- Note: at this point rho is constant, so it could be taken out
//        i.e. rho*( ddt(T) + div(phi, T) ).
//        However, in the event of adding and EOS to link rho, P and T,
//        the solver must become based on rhoPhi instead of phi,
//        so that it keeps its conservative characteristics.

fvScalarMatrix TEqn
(
    cp.F()*fvm::ddt(constEq.rho(), T)
  + cp.F()*fvm::div(fvc::interpolate(constEq.rho())*phi, T, "div(phi,T)")
  - fvm::laplacian(kt.F(), T, "laplacian(kt,T)")
  ==
    (
        newtonDissSwitch
        ?
        //constEq.tauN() && symm(fvc::grad(U))  // can be simplified...
        constEq.tauN() && fvc::grad(U)  // T:S==T:grad(U) for symmetric T
        :
        tmp<volScalarField> ( zeroDissTerm )
    )
  + (
        velastDissSwitch
        ?
        constEq.energyExtraTerms()
        :
        tmp<volScalarField> ( zeroDissTerm )
    )
  //- fvm::SuSp(-cte*fvc::DDt(p), T)  // eventually when an EOS is added
  + fvOptions(cp.F()*constEq.rho(), T)
);

TEqn.relax();

fvOptions.constrain(TEqn);

TEqn.solve();

fvOptions.correct(T);

#include "updateThermFuns.H"

