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
    rhoFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"

#include "fvMatrix.H"

#include "surfaceInterpolation.H"
#include "LimitedScheme.H"
#include "Minmod.H"

#include "fluxScheme.H"
#include "ePsiThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
	Info<< "\n---------Version 0.0--------\n" << endl;
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
	#include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //#include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
    	const float t0 = runTime.elapsedCpuTime();
    	Info << "start:" << runTime.cpuTimeIncrement() << endl;

    	//store previous iterations for relaxation
    	rho.storePrevIter();
    	rhoU.storePrevIter();
    	rhoE.storePrevIter();

    	volScalarField c(sqrt(thermo.Cp()/(thermo.Cv()*psi)));
    	surfaceScalarField phiv(fvc::interpolate(U) & mesh.Sf());
    	surfaceScalarField cSf(fvc::interpolate(c)*mesh.magSf());
    	surfaceScalarField amaxSf(max(mag(phiv + cSf),mag(phiv - cSf)));

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

    	Info << runTime.cpuTimeIncrement() << endl;
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << runTime.cpuTimeIncrement() << endl;
        tmp<fv::fluxScheme<scalar> > flux = fv::fluxScheme<scalar>::New(mesh,thermo,U,rho,mesh.divScheme("flux"));
        Info << "1-- " << runTime.cpuTimeIncrement() << endl;
        //flux().calculate(pAve,rhoAve,UAve);
        flux().flux(rhoFlux,UFlux,EFlux);
        Info << "2-- " << runTime.cpuTimeIncrement() << endl;
        Info << runTime.cpuTimeIncrement() << endl;

        volScalarField muEff(turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U)))); //part of equation.4

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(rhoFlux));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(UFlux));

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT());

        if (!inviscid)
        {
        	fvMatrix<vector> UEquation(
        		fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
        	);
        	UEquation.relax();
            solve(UEquation);   //equation.2


            // TODO verify if the bound function exists
            double Umax;
            double Umin;
            mesh.schemesDict().readIfPresent("Umax", Umax);
            mesh.schemesDict().readIfPresent("Umin", Umin);
            const vector UminVect(Umin,Umin,Umin);
            const vector UmaxVect(Umax,Umax,Umax);
			if((max(U).value()[0] > Umax) | (max(U).value()[1] > Umax))
			{
				Info << "U limited to " << Umax << " it was " << max(U).value() << endl;
				U.min(dimensionedVector("Umax", U.dimensions(), UmaxVect));
			}
			if((min(U).value()[0] < Umin) | (min(U).value()[1] < Umin))
			{
				Info << "U limited to " << Umin << " it was " << min(U).value() << endl;
				U.max(dimensionedVector("Umin", U.dimensions(), UminVect));
			}
            rhoU = rho*U;
        }

        // --- Solve energy
        sigmaDotU = (
                	fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              	  + (mesh.Sf() & fvc::interpolate(tauMC))
            		)
            		& fvc::interpolate(U); //UAve; //(a_pos*U_pos + a_neg*U_neg); TODO verify if this might be a problem for stability

        solve
		(
			fvm::ddt(rhoE)
			+ fvc::div(EFlux)
			- fvc::div(sigmaDotU)
		);

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();

		bool eLimited = false;
		scalarField Tmin(1);
		scalarField Tmax(1);
		labelList cell(1);
		cell[0] = 0;
		mesh.schemesDict().readIfPresent("Tmax", Tmax[0]);
		mesh.schemesDict().readIfPresent("Tmin", Tmin[0]);
		scalar emin  = thermo.e(Tmin,cell)()[0];
		scalar emax  = thermo.e(Tmax,cell)()[0];
		if(min(e).value()<emin)
		{
			Info << "e limited to " << emin << " it was " << min(e).value() << endl;
			e.max(emin);
			eLimited = true;
		}
		if(max(e).value()>emax)
		{
			Info << "e limited to " << emax << " it was " << max(e).value() << endl;
			e.min(emax);
			eLimited = true;
		}
		thermo.correct();

		if(eLimited)
		{
			rhoE = rho*(e + 0.5*magSqr(U));
		}
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        double rhomax;
        double rhomin;
        mesh.schemesDict().readIfPresent("rhomax", rhomax);
        mesh.schemesDict().readIfPresent("rhomin", rhomin);
        if (!inviscid)
        {
        	//volScalarField& e0 = e.oldTime();
        	//e0 = e;
            volScalarField k("k", thermo.Cp()*muEff/Pr);
            fvMatrix<double> eEquation(
            	fvm::ddt(rho, e) - fvc::ddt(rho, e)
			  - fvm::laplacian(turbulence->alphaEff()-turbulence->alpha(), e)
			  - fvm::laplacian(k/thermo.Cv(), e)
            );
            eEquation.relax();
			solve(eEquation);
            thermo.correct();

            if(min(rho).value()<rhomin)
			{
				Info << "rho limited to " << rhomin << " it was " << min(rho).value() << endl;
				rho.max(rhomin);
			}
            if(max(rho).value()>rhomax)
			{
				Info << "rho limited to " << rhomax << " it was " << max(rho).value() << endl;
				rho.min(rhomax);
			}
            //e.relax();
            rhoE = rho*(e + 0.5*magSqr(U));
        }

        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();

        rho.boundaryField() = psi.boundaryField()*p.boundaryField();
        if(min(rho).value()<rhomin)
		{
			Info << "rho limited to " << rhomin << " it was " << min(rho).value() << endl;
			rho.max(rhomin);
		}
		if(max(rho).value()>rhomax)
		{
			Info << "rho limited to " << rhomax << " it was " << max(rho).value() << endl;
			rho.min(rhomax);
		}

        turbulence->correct();

        runTime.write();

        Info << runTime.cpuTimeIncrement() << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "   Time/iter = " << runTime.elapsedCpuTime() - t0 << " s"
            << " iteration = " << runTime.timeIndex() << " "
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
