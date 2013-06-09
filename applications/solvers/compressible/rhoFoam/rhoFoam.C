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

        // --- upwind interpolation of primitive fields on faces
    	// TODO the kurganov flux splitting implementation shall be hard coded in the finiteVolume library
        surfaceInterpolation::debug = false;

    	surfaceScalarField rho_pos
        (
            fvc::extrapolate(rho, pos, "reconstruct(rho)")
        );
    	surfaceInterpolation::debug = false;

        surfaceScalarField rho_neg
        (
            fvc::extrapolate(rho, neg, "reconstruct(rho)")
        );


        surfaceVectorField rhoU_pos
        (
            fvc::extrapolate(rhoU, pos, "reconstruct(U)")
        );

        surfaceVectorField rhoU_neg
        (
            fvc::extrapolate(rhoU, neg, "reconstruct(U)")
        );

        volScalarField rPsi(1.0/psi);

        surfaceScalarField rPsi_pos
        (
            fvc::extrapolate(rPsi, pos, "reconstruct(T)")
        );

        surfaceScalarField rPsi_neg
        (
            fvc::extrapolate(rPsi, neg, "reconstruct(T)")
        );


        /*surfaceScalarField e_pos
        (
            fvc::extrapolate(e, pos, "reconstruct(T)")
        );

        surfaceScalarField e_neg
        (
            fvc::extrapolate(e, neg, "reconstruct(T)")
        );*/

        surfaceVectorField U_pos(rhoU_pos/rho_pos);
        surfaceVectorField U_neg(rhoU_neg/rho_neg);

        surfaceScalarField p_pos(rho_pos*rPsi_pos);  //rPsi=p/rho ==> rPsi=RT (perfect gas assumption)
        surfaceScalarField p_neg(rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos(U_pos & mesh.Sf());
        surfaceScalarField phiv_neg(U_neg & mesh.Sf());

        volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));  //perfect gas assumption

        surfaceScalarField cSf_pos
        (
            fvc::extrapolate(c, pos, "reconstruct(T)")*mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            fvc::extrapolate(c, neg, "reconstruct(T)")*mesh.magSf()
        );

        surfaceScalarField ap
        ("ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)  //equation.8.a
        );
        surfaceScalarField am
        ("am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)  //(-1)*equation.8.b
        );

        surfaceScalarField a_pos(ap/(ap - am));  //equation.9.b (i.e. Kurganov)

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf(am*a_pos);  //(-1)*equation.10.b (i.e. Kurganov)

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;  //(-1)*equation.10.a (i.e. Tadmor)
            a_pos = 0.5;        //equation.9.a (i.e. Tadmor)
        }

        surfaceScalarField a_neg(1.0 - a_pos);

        surfaceScalarField phiv_pos1("phiv_pos",phiv_pos);
        surfaceScalarField phiv_neg1("phiv_neg",phiv_neg);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos",phiv_pos - aSf);
        surfaceScalarField aphiv_neg(phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg; //equation.7 for rho

        //volScalarField phi_(fvc::div(phiv_pos1,rho));
        //Info << max(mag(fvc::div(phi) - phi_)).value() << endl;

        Istream& Is(mesh.divScheme("kurg"));
        word dummy(Is);
        fv::kurganovConvectionScheme<double> kurg(mesh,pos,Is);
        Istream& Is1(mesh.divScheme("kurgV"));
        word dummy1(Is1);
        fv::kurganovConvectionScheme<vector> kurgV(mesh,pos,Is1);
        //volScalarField phi_(kurg.fvcDiv(phiv_pos1,phiv_neg1,rho));

        fv::goudonovFluxScheme sTest(mesh,phi,p,U,rho);
        sTest.calculate(pAve,rhoAve,UAve);



        surfaceScalarField rhoFlux(rhoAve*(UAve & mesh.Sf()));
        surfaceVectorField UFlux(rhoFlux*UAve + pAve*mesh.Sf());
        surfaceScalarField EFlux((UAve & mesh.Sf())*(pAve/(1.4-1) + 0.5*rhoAve*magSqr(UAve) + pAve));

		//#include "cellDebug.H"

        volScalarField muEff(turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U)))); //part of equation.4

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(rhoFlux));
        //solve(fvm::ddt(rho) + fvc::div(phiv_pos1,rho)); //equation.1
        //solve(fvm::ddt(rho) + kurg.fvcDiv(phiv_pos1,phiv_neg1,rho)); //equation.1
        //solve(fvm::ddt(rho) + fvc::div(phi)); //equation.1
        //solve(fvm::ddt(rho) + fvc::div(phi1) + divphi2); //equation.1
        rho.relax();

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(UFlux));
        //solve(fvm::ddt(rhoU) + fvc::div(phiv_pos1,rhoU)
        //        		+ fvc::div(kurg.interpolate(a_pos,a_neg,p)*mesh.Sf())); //equation.2 inviscid
        /*solve(fvm::ddt(rhoU) + kurgV.fvcDiv(phiv_pos1,phiv_neg1,rhoU)
        		+ fvc::div(kurg.interpolate(a_pos,a_neg,p)*mesh.Sf()));*/ //equation.2 inviscid
        //solve(fvm::ddt(rhoU) + fvc::div(phiUp)); //equation.2 inviscid
        //solve(fvm::ddt(rhoU) + fvc::div(phiUp1) + divphiUp2); //equation.2 inviscid
        rhoU.relax();

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT());

        if (!inviscid)
        {
        	//volVectorField& U0 = U.oldTime();
        	//U0 = U;
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
            		& (a_pos*U_pos + a_neg*U_neg);

        solve
		(
			fvm::ddt(rhoE)
			+ fvc::div(EFlux)
			- fvc::div(sigmaDotU)
		);
        /*solve
		(
			fvm::ddt(rhoE)
			+ fvc::div(phiv_pos1,rho*(e + 0.5*magSqr(U)) + p,"kurg")
    		+ fvc::div(kurg.interpolate(aSf,-aSf,p))
			- fvc::div(sigmaDotU)
		);*/
        /*solve
		(
			fvm::ddt(rhoE)
			+ fvc::div(phiEp)
			- fvc::div(sigmaDotU)
		);*/
        /*solve
        (
			fvm::ddt(rhoE)
			+ fvc::div(phiEp1)
			+ divphiEp2
			- fvc::div(sigmaDotU)
        );*/
        rhoE.relax();

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

        /*phiEp_e = aphiv_pos*rho_pos*e_pos + aphiv_neg*rho_neg*e_neg;
        phiEp_U = aphiv_pos*rho_pos*0.5*magSqr(U_pos) + aphiv_neg*rho_neg*0.5*magSqr(U_neg);
        phiEp_p = aphiv_pos*p_pos + aphiv_neg*p_neg + aSf*p_pos - aSf*p_neg;
        divphiEp = fvc::div(phiEp);
        divphi = fvc::div(phi);
        diffT = (-fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(e) +
        	fvc::interpolate(turbulence->alpha())*fvc::snGrad(e)
        	- fvc::interpolate(thermo.Cp()*muEff/Pr)*fvc::snGrad(T))*mesh.magSf();
        ddtrho = fvc::ddt(rho);
        ddtrhoE = fvc::ddt(rho,e);
        divphiEp1 = fvc::div(phiEp);*/

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
