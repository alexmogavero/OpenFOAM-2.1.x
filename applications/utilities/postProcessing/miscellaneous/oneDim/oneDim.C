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
#include "sampledPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject Theader
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject pheader
    (
    	"p",
    	runTime.timeName(),
    	mesh,
    	IOobject::MUST_READ
    );

    volScalarField T(Theader,mesh);
    volScalarField p(pheader,mesh);
    volVectorField U(Uheader,mesh);

    int N = 100;
    double Xmin = 0;
    double Xmax = 8;
    double DX = (Xmax - Xmin)/(N-1);

    scalarList X1dim;
    scalarList T1dim;
    scalarList p1dim;
    scalarList U1dim;

    for(int j=0; j<N; j++){
    	double X = Xmin + j*DX;

    	sampledPlane pl
    	(
    			"plane",
    			mesh,
    			plane(vector(X,0,0),vector(1,0,0)),
    			word::null
    	);

    	pl.update();
    	if(pl.area()>0){
    		tmp<scalarField> Ts = pl.sample(T);
    		tmp<scalarField> ps = pl.sample(p);
    		tmp<vectorField> Us = pl.sample(U);
    		const scalarField& Area = pl.magSf();


    		scalar Tave = 0;
    		scalar pave = 0;
    		scalar Uave = 0;
    		forAll(Ts(),i){

    			Tave = Tave + Area[i]*Ts()[i];
    			pave = pave + Area[i]*ps()[i];
    			Uave = Uave + Area[i]*Us()[i].x();
    		}

    		Tave = Tave/pl.area();
    		pave = pave/pl.area();
    		Uave = Uave/pl.area();

    		X1dim.append(X);
    		T1dim.append(Tave);
    		p1dim.append(pave);
    		U1dim.append(Uave);
    	}
    }

    if(writeResults){
    	Info << "Writing output to file" << endl;

    	OFstream outFile("oneDim.dat");
    	outFile << "X p T U_x" << endl;
    	forAll(X1dim,i){
    		outFile << X1dim[i] << " " << p1dim[i] << " " << T1dim[i] << " " << U1dim[i] <<  endl;
    	}
    }else{
    	forAll(X1dim,i){
    		Info << X1dim[i] << " " << p1dim[i] << " " << T1dim[i] << " " << U1dim[i] <<  endl;
    	}
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
