/*
 * debugFields.H
 *
 *  Created on: 31 Jan 2013
 *      Author: trb12187
 */

surfaceScalarField phiEp_e
(
	IOobject
	(
		"phiEp_e",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("phiEp_e", cmptMultiply(phi.dimensions(),pow(U.dimensions(),2)), 0.0)
);
surfaceScalarField phiEp_U
(
	IOobject
	(
		"phiEp_U",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("phiEp_U", phiEp_e.dimensions(), 0.0)
);
surfaceScalarField phiEp_p
(
	IOobject
	(
		"phiEp_p",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("phiEp_p", phiEp_e.dimensions(), 0.0)
);
surfaceScalarField sigmaDotU
(
	IOobject
	(
		"sigmaDotU",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("sigmaDotU", phiEp_e.dimensions(), 0.0)
);
volScalarField divphiEp
(
	IOobject
	(
		"divphiEp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("divphiEp", cmptDivide(phiEp_e.dimensions(),mesh.V().dimensions()), 0.0)
);
volScalarField divphiEp1
(
	IOobject
	(
		"divphiEp1",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("divphiEp1", divphiEp.dimensions(), 0.0)
);
volScalarField divphiEp2
(
	IOobject
	(
		"divphiEp2",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("divphiEp2", divphiEp.dimensions(), 0.0)
);
volScalarField divphi
(
	IOobject
	(
		"divphi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("divphi", cmptDivide(phi.dimensions(),mesh.V().dimensions()), 0.0)
);
surfaceScalarField diffT
(
	IOobject
	(
		"diffT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("diffT", phiEp_e.dimensions(), 0.0)
);
volScalarField ddtrho
(
	IOobject
	(
		"ddtrho",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("ddtrho", divphi.dimensions(), 0.0)
);
volScalarField ddtrhoE
(
	IOobject
	(
		"ddtrhoE",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("ddtrhoE", divphiEp.dimensions(), 0.0)
);





