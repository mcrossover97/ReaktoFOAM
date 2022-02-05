Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar eps0
(
	transportProperties.lookup("eps0")
);

dimensionedScalar k0
(
	transportProperties.lookup("k0")
);

dimensionedScalar mu
(
	transportProperties.lookup("mu")
);

dimensionedScalar temp
(
	transportProperties.lookup("temp")
);

dimensionedScalar dispCo
(
	transportProperties.lookup("dispCo")
);

dimensionedScalar xNo
(
	transportProperties.lookup("xNo")
);

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField De
(
    IOobject
    (
        "De",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	dimensionedTensor("De",dimensionSet(0,2,-1,0,0,0,0),tensor(dispCo.value(),0,0,0,dispCo.value(),0,0,0,dispCo.value()))
);

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

#include "createPhi.H"
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());