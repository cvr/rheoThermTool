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


//fvSchemes (static_cast<const objectRegistry&>(*this))


//Foam::thermFunModels thermFunModels
//(
//    Foam::IOobject
//    (
//        Foam::fvMesh::defaultRegion,
//        runTime.timeName(),
//        runTime,
//        Foam::IOobject::MUST_READ
//    )
//);


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

/*
thermFunModel mu("mu", T);
thermFunModel ert("ert", T);

thermFunModel cp("cp", T);
thermFunModel k("k", T);
*/




//tensor gradU_(cttProperties.subDict("rheoTestFoamParameters").lookup("gradU"));
//scalarField gammaEpsilonDotL_(cttProperties.subDict("rheoTestFoamParameters").lookup("gammaEpsilonDotL"));

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


Info << "»»»»»»»»» DEBUG «««««««««« a0" << endl;

// Create constitutive equation
constitutiveModel constEq(U, phi, T);

Info << "»»»»»»»»» DEBUG «««««««««« a1" << endl;

bool isGNF(constEq.isGNF());

word CM = cttProperties.subDict("parameters").lookup("type");
scalar lambdaMax(0.0);
scalar lambdaMin(1e20);

if (!isGNF) // doesn't care about ramp, although lambdaMax/Min have no sense if ramp is false
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


