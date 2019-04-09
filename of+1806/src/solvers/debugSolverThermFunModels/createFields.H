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


Info<< "Creating thermal functions\n" << endl;
thermFunModel cp("cp", T);
thermFunModel kt("kt", T);
thermFunModel eta("eta", T);
thermFunModel lambda("lambda", T);

