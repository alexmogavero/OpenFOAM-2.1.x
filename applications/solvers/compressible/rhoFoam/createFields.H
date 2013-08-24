Info<< "Reading thermophysical properties\n" << endl;

autoPtr<basicPsiThermo> pThermo
(
    basicPsiThermo::New(mesh)
);
basicPsiThermo& thermo = pThermo();

volScalarField& p = thermo.p();
volScalarField& e = thermo.e();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.internalField()) > 0.0)
{
    inviscid = false;
}

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

#include "rhoBoundaryTypes.H"
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho(),
    rhoBoundaryTypes
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);


surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

/*surfaceScalarField pAve
(
	IOobject
	(
		"pAve",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("pAve", p.dimensions(), 0.0)
);

surfaceScalarField rhoAve
(
	IOobject
	(
		"rhoAve",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("rhoAve", rho.dimensions(), 0.0)
);

surfaceVectorField UAve
(
	IOobject
	(
		"UAve",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("UAve", U.dimensions(), vector::zero)
);

surfaceScalarField eAve
(
	IOobject
	(
		"eAve",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("eAve", e.dimensions(), 0.0)
);

surfaceScalarField TAve
(
	IOobject
	(
		"TAve",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("TAve", T.dimensions(), 0.0)
);*/

surfaceScalarField rhoFlux
(
	IOobject
	(
		"rhoFlux",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("rhoFlux", rho.dimensions()*U.dimensions()*mesh.Sf().dimensions(), 0.0)
);

surfaceScalarField EFlux
(
	IOobject
	(
		"EFlux",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("EFlux", rhoFlux.dimensions()*e.dimensions(), 0.0)
);

surfaceVectorField UFlux
(
	IOobject
	(
		"UFlux",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedVector("UFlux", rhoFlux.dimensions()*U.dimensions(), vector::zero)
);


#include "debugFields.H"

