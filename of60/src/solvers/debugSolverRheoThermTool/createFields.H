Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Preparing IO for constant/constitutiveProperties\n" << endl;
IOdictionary cttProperties
(
    IOobject
    (
        "constitutiveProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    )
);


Info<< "Reading field T\n" << endl;
dimensionedScalar T0 = cttProperties.subDict("thermophysicalProperties").lookup("T0");
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("T", T0.dimensions(), T0.value())
);


Info<< "Checking existance of passive scalar\n" << endl;
bool sPS = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false);

autoPtr<volScalarField> C
(
    sPS != true
    ?
    NULL
    :
    new volScalarField
    (
	IOobject
	(
	    "C",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    )
);


Info<< "Creating model(s) for constitutive equation\n" << endl;
constitutiveModel constEq(U, phi, T);


Info<< "Creating thermal functions\n" << endl;

Info<< "Setting thermFun for cp (isobaric specific heat)" << endl;
thermFunModel cp("cp", T);
cp.calcThermFun(T, T0);

Info<< "Setting thermFun for kt (thermal conductivity)" << endl;
thermFunModel kt("kt", T);
kt.calcThermFun(T, T0);

//- Comment: dirty way to do it but if I create
//           thermFunModel k("k", T) inside the
//           if statement, it gets within that
//           scope only and I cannot get it out!
//           This works as a kind of hack.

thermFunModel* k = nullptr;
if (&constEq.kRef() != nullptr)
{
    Info<< "Setting thermFun for k (GNF)" << endl;
    k = new thermFunModel("k", T);
    (*k).calcThermFun(T, T0);
}

thermFunModel* etaMin = nullptr;
if (&constEq.etaMinRef() != nullptr)
{
    Info<< "Setting thermFun for etaMin (GNF)" << endl;
    etaMin = new thermFunModel("etaMin", T);
    (*etaMin).calcThermFun(T, T0);
}

thermFunModel* etaMax = nullptr;
if (&constEq.etaMaxRef() != nullptr)
{
    Info<< "Setting thermFun for etaMax (GNF)" << endl;
    etaMax = new thermFunModel("etaMax", T);
    (*etaMax).calcThermFun(T, T0);
}

thermFunModel* etaS = nullptr;
if (&constEq.etaSRef() != nullptr)
{
    Info<< "Setting thermFun for etaS (VE)" << endl;
    etaS = new thermFunModel("etaS", T);
    (*etaS).calcThermFun(T, T0);
}

thermFunModel* etaP = nullptr;
if (&constEq.etaPRef() != nullptr)
{
    Info<< "Setting thermFun for etaP (VE)" << endl;
    etaP = new thermFunModel("etaP", T);
    (*etaP).calcThermFun(T, T0);
}

thermFunModel* lambda = nullptr;
if (&constEq.lambdaRef() != nullptr) {
    Info<< "Setting thermFun for lambda (VE)" << endl;
    lambda = new thermFunModel("lambda", T);
    (*lambda).calcThermFun(T, T0);
}




Info << "»»»»»»»»» DEBUG «««««««««« a1" << endl;

word CM = cttProperties.subDict("parameters").lookup("type");
scalar lambdaMax(0.0);
scalar lambdaMin(1e20);

if (!constEq.isGNF()) // doesn't care about ramp, although lambdaMax/Min have no sense if ramp is false
{
    if(CM=="multiMode")
    {
        PtrList<entry> modelEntries(cttProperties.subDict("parameters").lookup("models"));
        forAll(modelEntries, modelI)
        {
            dimensionedScalar lambdaI(modelEntries[modelI].dict().lookup("lambda"));
            lambdaMax = Foam::max(lambdaMax, lambdaI.value());
            lambdaMin = Foam::min(lambdaMin, lambdaI.value());    
        }
    }
    else
    {
        dimensionedScalar lambdatmp(cttProperties.subDict("parameters").lookup("lambda"));
        lambdaMax = lambdatmp.value();
        lambdaMin = lambdaMax;
    }
}



Info << "»»»»»»»»» DEBUG «««««««««« a2" << endl;

//volVectorField centres = mesh.C();
//volScalarField volume= Sj.mesh().V().field();
volScalarField X = mesh.C().component(vector::X);
volScalarField Y = mesh.C().component(vector::Y);
volScalarField Xadim = sqrt(sqr(X / mag(max(X) - min(X))));
volScalarField Yadim = sqrt(sqr(Y / mag(max(Y) - min(Y))));
volScalarField XYadim = sqrt(sqr(Xadim) + sqr(Yadim));

Info << "»»»»»»»»» DEBUG «««««««««« a3" << endl;



