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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


#include "createMRF.H"
#include "createFvOptions.H"


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


#include "createThermFuns.H"


