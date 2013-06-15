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
    rhoCentralFoam

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
#include "kurganovConvectionScheme.H"

#include "exactReinmannSolver.H"
#include "wafFluxScheme.H"
#include "goudonovFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
const Foam::word fv::exactReinmannSolver::typeName = "ExactReinmann";
const Foam::word fv::wafFluxScheme::typeName = "WAF";
const Foam::word fv::goudonovFluxScheme::typeName = "WAF";
int main(int argc, char *argv[])
{
	Info<< "\n---------Version 2.0--------\n" << endl;
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
	#include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
    	const float t0 = runTime.elapsedCpuTime();
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

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        fv::goudonovFluxScheme flux(mesh,phi,p,U,rho);
        flux.calculate(pAve,rhoAve,UAve);

        surfaceScalarField rhoFlux(rhoAve*(UAve & mesh.Sf()));
        surfaceVectorField UFlux(rhoFlux*UAve + pAve*mesh.Sf());
        dimensionedScalar e0("e0",e.dimensions(),thermo.e(0*T,1)()[0]);
        surfaceScalarField EFlux((UAve & mesh.Sf())*(e0*rhoAve + pAve/(1.4-1) + 0.5*rhoAve*magSqr(UAve) + pAve));

		#include "cellDebug.H"

        volScalarField muEff(turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U)))); //part of equation.4

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(rhoFlux));
        //rho.relax();

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(UFlux));
        //rhoU.relax();

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
            		& UAve; //(a_pos*U_pos + a_neg*U_neg);

        solve
		(
			fvm::ddt(rhoE)
			+ fvc::div(EFlux)
			- fvc::div(sigmaDotU)
		);
        //rhoE.relax();

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();

        //TODO limit the energy starting from the min and max Temperature
        double emax;
		double emin;
		mesh.schemesDict().readIfPresent("emax", emax);
		mesh.schemesDict().readIfPresent("emin", emin);
		if(min(e).value()<emin)
		{
			Info << "e limited to " << emin << " it was " << min(e).value() << endl;
			e.max(emin);
		}
		if(max(e).value()>emax)
		{
			Info << "e limited to " << emax << " it was " << max(e).value() << endl;
			e.min(emax);
		}
		thermo.correct();

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

        //Info << p[0] << " ";
        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        //Info << p[0] << " " << (e[0]*(1.4-1)) << endl;

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

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "   Time/iter = " << runTime.elapsedCpuTime() - t0 << " s"
            << " iteration = " << runTime.timeIndex() << " "
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
