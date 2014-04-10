/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Mach

Description
    Calculates and optionally writes the local Mach number from the velocity
    field U at each time.

    The -nowrite option just outputs the max value without writing the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "basicThermo.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject UAveHeader
	(
		"UAve",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ
	);

	IOobject rhoAveHeader
	(
		"rhoAve",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ
	);

    // Check U and rho exists
    if (UHeader.headerOk() && rhoHeader.headerOk())
    {
        volVectorField U(UHeader, mesh);
        volScalarField rho(rhoHeader, mesh);

        surfaceScalarField mfr
        (
        	IOobject
            (
                "mfr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (mesh.Sf() & fvc::interpolate(U))*fvc::interpolate(rho)
        );

        surfaceScalarField mfrAve
		(
			IOobject
			(
				"mfrAve",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("mfrAve", mfr.dimensions(), 0.0)
		);

        if (UAveHeader.headerOk() && rhoAveHeader.headerOk())
        {
        	surfaceVectorField UAve(UAveHeader, mesh);
        	surfaceScalarField rhoAve(rhoAveHeader, mesh);
        	mfrAve = (mesh.Sf() & UAve)*rhoAve;
        }

        scalar mfrTot = 0;
        scalar mfrAveTot = 0;

        Info<< "\nMass Flow rate [kg/s]" << endl;
        forAll(mfr.boundaryField(),patchi)
        {
        	scalar mfrPatch = gSum(mfr.boundaryField()[patchi]);
        	mfrTot = mfrTot + mfrPatch;
        	Info<< mesh.boundary()[patchi].name() << " " << mfrPatch;

        	if (UAveHeader.headerOk() && rhoAveHeader.headerOk())
        	{
        		scalar mfrAvePatch = gSum(mfrAve.boundaryField()[patchi]);
        		mfrAveTot = mfrAveTot + mfrAvePatch;
        		Info << " " << mfrAvePatch;
        	}

			Info << endl;
        }

        Info << "-----------------------------------------------" << endl;
        Info << "TOT " << mfrTot;
        if (UAveHeader.headerOk() && rhoAveHeader.headerOk())
        {
        	Info << " " << mfrAveTot;
        }
        Info << endl;

        if (writeResults)
        {
            mfr.write();

            if (UAveHeader.headerOk() && rhoAveHeader.headerOk())
			{
            	mfrAve.write();
			}
        }
    }
    else
    {
        Info<< "    Missing U or rho" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
