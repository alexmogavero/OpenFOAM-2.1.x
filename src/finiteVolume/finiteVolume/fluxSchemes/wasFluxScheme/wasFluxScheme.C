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

\*---------------------------------------------------------------------------*/

#include "wasFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void
wasFluxScheme<Type>::calculate
(
	surfaceScalarField& pAve,
	surfaceScalarField& rhoAve,
	surfaceVectorField& UAve
)
{
	riemann();

	const fvMesh& mesh = p_.mesh();
	const labelUList& owner = mesh.owner();
	const labelUList& neighbour = mesh.neighbour();
	const volVectorField& C = mesh.C();
	const volVectorField gradRho0(fvc::grad(rho0_));
	const volVectorField gradRho1(fvc::grad(rho1_));
	const volVectorField gradRho2(fvc::grad(rho2_));
	const volVectorField gradRho3(fvc::grad(rho3_));
	const volTensorField gradU(fvc::grad(U_));

	forAll(p1_,face)
	{
		label own = owner[face];
		label nei = neighbour[face];
		const vector d = C[nei] - C[own];
		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;
		const vector vL = U_[own] - uL*n;
		const vector vR = U_[nei] - uR*n;
		const tensor gradVL = gradU[own] - (gradU[own] & n)*n;
		const tensor gradVR = gradU[nei] - (gradU[nei] & n)*n;

		double w[3];
		const double rho[4] = {rho0_[face], rho1_[face], rho2_[face], rho3_[face]};
		const double speed[3] = {speed0_[face], speed1_[face], speed2_[face]};

		const double rhoL[4] = {rho0_[face]-(d & gradRho0[own]), rho1_[face]-(d & gradRho1[own]),
				rho2_[face]-(d & gradRho2[own]), rho3_[face]-(d & gradRho3[own])};
		const double rhoR[4] = {rho0_[face]+(d & gradRho0[nei]), rho1_[face]+(d & gradRho1[nei]),
				rho2_[face]+(d & gradRho2[nei]), rho3_[face]+(d & gradRho3[nei])};

		weight(w, rho, speed, rhoL, rhoR, d);

		const double p[4] = {p0_[face], p1_[face], p2_[face], p3_[face]};
		pAve[face] = average(p,w);
		rhoAve[face] = average(rho,w);
		const double u[4] = {u0_[face], u1_[face], u2_[face], u3_[face]};
		const scalar uAve = average(u,w);

		vector vAve = averageV(speed1_[face],vL,vR,gradVL,gradVR,d);
		UAve[face] = uAve*n + vAve;

	}

	const surfaceScalarField::GeometricBoundaryField& bP = p1_.boundaryField();
	forAll(bP, patchi)
	{
		const scalarField p0 = p0_.boundaryField()[patchi];
		const scalarField p1 = p1_.boundaryField()[patchi];
		const scalarField p2 = p2_.boundaryField()[patchi];
		const scalarField p3 = p3_.boundaryField()[patchi];
		const scalarField rho0 = rho0_.boundaryField()[patchi];
		const scalarField rho0B = rho0B_[patchi];
		const scalarField rho1 = rho1_.boundaryField()[patchi];
		const scalarField rho1B = rho1B_[patchi];
		const scalarField rho2 = rho2_.boundaryField()[patchi];
		const scalarField rho2B = rho2B_[patchi];
		const scalarField rho3 = rho3_.boundaryField()[patchi];
		const scalarField rho3B = rho3B_[patchi];
		const vectorField UL = U_.boundaryField()[patchi].patchInternalField();
		const vectorField UR = 2*U_.boundaryField()[patchi] - UL;
		const scalarField u0 = u0_.boundaryField()[patchi];
		const scalarField u1 = u1_.boundaryField()[patchi];
		const scalarField u2 = u2_.boundaryField()[patchi];
		const scalarField u3 = u3_.boundaryField()[patchi];
		const scalarField speed0 = speed0_.boundaryField()[patchi];
		const scalarField speed1 = speed1_.boundaryField()[patchi];
		const scalarField speed2 = speed2_.boundaryField()[patchi];

		const tensorField gradUL = gradU.boundaryField()[patchi].patchInternalField();
		const tensorField gradUR = 2*gradU.boundaryField()[patchi] - gradU.boundaryField()[patchi].patchInternalField();
		const vectorField gradRho0L = gradRho0.boundaryField()[patchi].patchInternalField();
		const vectorField gradRho1L = gradRho1.boundaryField()[patchi].patchInternalField();
//		const vectorField gradRho1R = 2*gradRho1.boundaryField()[patchi] - gradRho1L;
		const vectorField gradRho2L = gradRho2.boundaryField()[patchi].patchInternalField();
//		const vectorField gradRho2R = 2*gradRho2.boundaryField()[patchi] - gradRho2L;
		const vectorField gradRho3L = gradRho3.boundaryField()[patchi].patchInternalField();

		const vectorField d = 2*mesh.boundary()[patchi].delta();

		scalarField& pAveB = pAve.boundaryField()[patchi];
		scalarField& rhoAveB = rhoAve.boundaryField()[patchi];
		vectorField& UAveB = UAve.boundaryField()[patchi];

		forAll(pAveB,face)
		{
			double w[3];
			const double pB[4] = {p0[face], p1[face], p2[face], p3[face]};
			const double rhoB[4] = {rho0[face], rho1[face], rho2[face], rho3[face]};
			const double uB[4] = {u0[face], u1[face], u2[face], u3[face]};
			const double speedB[3] = {speed0[face], speed1[face], speed2[face]};

			const double rhoLB[4] = {rho0[face]-(d[face] & gradRho0L[face]), rho1[face]-(d[face] & gradRho1L[face]),
					rho2[face]-(d[face] & gradRho2L[face]), rho3[face]-(d[face] & gradRho3L[face])};
			const double rhoRB[4] = {rho0B[face],rho1B[face],rho2B[face],rho3B[face]};

			const vector n = mesh.boundary()[patchi].Sf()[face]/mesh.boundary()[patchi].magSf()[face];

			weight(w, rhoB, speedB, rhoLB, rhoRB, d[face]);


			pAveB[face] = average(pB,w);
			rhoAveB[face] = average(rhoB,w);
			const scalar uAve = average(uB,w);

			const vector vL = UL[face] - (UL[face] & n)*n;
			const vector vR = UR[face] - (UR[face] & n)*n;
			const tensor gradVL = gradUL[face] - (gradUL[face] & n)*n;
			const tensor gradVR = gradUR[face] - (gradUR[face] & n)*n;
			vector vAve = averageV(speed1[face],vL,vR,gradVL,gradVR,d[face]);
			UAveB[face] = uAve*n + vAve;

		}
	}
}

template<class Type>
void
wasFluxScheme<Type>::riemann
(
)
{
	const fvMesh& mesh = p_.mesh();
	const labelUList& owner = mesh.owner();
	const labelUList& neighbour = mesh.neighbour();
	const tmp<volScalarField> Cp = thermo_.Cp();
	const tmp<volScalarField> Cv = thermo_.Cv();
	const surfaceScalarField gamma(fvc::interpolate(Cp()/Cv()));

	forAll(p1_,face)
	{
		label own = owner[face];
		label nei = neighbour[face];

		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;

		exactRiemannSolver riemann(p_[own],p_[nei],rho_[own],rho_[nei],uL,uR,gamma[face]);
		riemann.solve();

		speed0_[face] = riemann.speed(0);
		speed1_[face] = riemann.speed(2);
		speed2_[face] = riemann.speed(4);
		double p[4];
		sampling(&exactRiemannSolver::p,riemann,p);
		p0_[face] = p[0];
		p1_[face] = p[1];
		p2_[face] = p[2];
		p3_[face] = p[3];
		double rho[4];
		sampling(&exactRiemannSolver::rho,riemann,rho);
		rho0_[face] = rho[0];
		rho1_[face] = rho[1];
		rho2_[face] = rho[2];
		rho3_[face] = rho[3];
		double u[4];
		sampling(&exactRiemannSolver::u,riemann,u);
		u0_[face] = u[0];
		u1_[face] = u[1];
		u2_[face] = u[2];
		u3_[face] = u[3];
	}

	const surfaceScalarField::GeometricBoundaryField& bP = p1_.boundaryField();
	setBoundary();
	forAll(bP, patchi)
	{
		const scalarField pL = p_.boundaryField()[patchi].patchInternalField();
		const scalarField pR = pB_[patchi];
		const scalarField pRR = pBB_[patchi];
		const scalarField rhoL = rho_.boundaryField()[patchi].patchInternalField();
		const scalarField rhoR = rhoB_[patchi];
		const scalarField rhoRR = rhoBB_[patchi];
		const vectorField UL = U_.boundaryField()[patchi].patchInternalField();
		const vectorField UR = UB_[patchi];
		const vectorField URR = UBB_[patchi];
		const scalarField gammab = gamma.boundaryField()[patchi];
		scalarField& p0B = p0_.boundaryField()[patchi];
		scalarField& p0BR = p0B_[patchi];
		scalarField& p1B = p1_.boundaryField()[patchi];
		scalarField& p1BR = p1B_[patchi];
		scalarField& p2B = p2_.boundaryField()[patchi];
		scalarField& p2BR = p2B_[patchi];
		scalarField& p3B = p3_.boundaryField()[patchi];
		scalarField& p3BR = p3B_[patchi];
		scalarField& rho0B = rho0_.boundaryField()[patchi];
		scalarField& rho0BR = rho0B_[patchi];
		scalarField& rho1B = rho1_.boundaryField()[patchi];
		scalarField& rho1BR = rho1B_[patchi];
		scalarField& rho2B = rho2_.boundaryField()[patchi];
		scalarField& rho2BR = rho2B_[patchi];
		scalarField& rho3B = rho3_.boundaryField()[patchi];
		scalarField& rho3BR = rho3B_[patchi];
		scalarField& u0B = u0_.boundaryField()[patchi];
		scalarField& u0BR = u0B_[patchi];
		scalarField& u1B = u1_.boundaryField()[patchi];
		scalarField& u1BR = u1B_[patchi];
		scalarField& u2B = u2_.boundaryField()[patchi];
		scalarField& u2BR = u2B_[patchi];
		scalarField& u3B = u3_.boundaryField()[patchi];
		scalarField& u3BR = u3B_[patchi];
		scalarField& speed0B = speed0_.boundaryField()[patchi];
		scalarField& speed0BR = speed0B_[patchi];
		scalarField& speed1B = speed1_.boundaryField()[patchi];
		scalarField& speed1BR = speed1B_[patchi];
		scalarField& speed2B = speed2_.boundaryField()[patchi];
		scalarField& speed2BR = speed2B_[patchi];

		forAll(p1B,face)
		{
			const vector n = mesh.boundary()[patchi].Sf()[face]/mesh.boundary()[patchi].magSf()[face];
			const double uL = UL[face] & n;
			const double uR = UR[face] & n;

			exactRiemannSolver riemann(pL[face],pR[face],rhoL[face],rhoR[face],uL,uR,gammab[face]);
			riemann.solve();

			speed0B[face] = riemann.speed(0);
			speed1B[face] = riemann.speed(2);
			speed2B[face] = riemann.speed(4);
			double p[4];
			sampling(&exactRiemannSolver::p,riemann,p);
			p0B[face] = p[0];
			p1B[face] = p[1];
			p2B[face] = p[2];
			p3B[face] = p[3];
			double rho[4];
			sampling(&exactRiemannSolver::rho,riemann,rho);
			rho0B[face] = rho[0];
			rho1B[face] = rho[1];
			rho2B[face] = rho[2];
			rho3B[face] = rho[3];
			double u[4];
			sampling(&exactRiemannSolver::u,riemann,u);
			u0B[face] = u[0];
			u1B[face] = u[1];
			u2B[face] = u[2];
			u3B[face] = u[3];


			const double uRR = URR[face] & n;

			exactRiemannSolver riemannR(pR[face],pRR[face],rhoR[face],rhoRR[face],uR,uRR,gammab[face]);
			riemannR.solve();

			speed0BR[face] = riemannR.speed(0);
			speed1BR[face] = riemannR.speed(2);
			speed2BR[face] = riemannR.speed(4);
			sampling(&exactRiemannSolver::p,riemannR,p);
			p0BR[face] = p[0];
			p1BR[face] = p[1];
			p2BR[face] = p[2];
			p3BR[face] = p[3];
			sampling(&exactRiemannSolver::rho,riemannR,rho);
			rho0BR[face] = rho[0];
			rho1BR[face] = rho[1];
			rho2BR[face] = rho[2];
			rho3BR[face] = rho[3];
			sampling(&exactRiemannSolver::u,riemannR,u);
			u0BR[face] = u[0];
			u1BR[face] = u[1];
			u2BR[face] = u[2];
			u3BR[face] = u[3];
		}
	}
}

template<class Type>
void
wasFluxScheme<Type>::setBoundary
(
)
{
	const fvMesh& mesh = p_.mesh();

	forAll(pB_, patchi)
	{
		scalarField& pB = pB_[patchi];
		scalarField& pBB = pBB_[patchi];
		vectorField& UB = UB_[patchi];
		vectorField& UBB = UBB_[patchi];
		scalarField& rhoB = rhoB_[patchi];
		scalarField& rhoBB = rhoBB_[patchi];

		if (isA<wallFvPatch>(mesh.boundary()[patchi]))
		{
			pB = p_.boundaryField()[patchi];
			pBB = p_.boundaryField()[patchi];
			rhoB = rho_.boundaryField()[patchi].patchInternalField();
			rhoBB = rho_.boundaryField()[patchi].patchInternalField();
			UB = 2*U_.boundaryField()[patchi] - U_.boundaryField()[patchi].patchInternalField();
			UBB = 3*UB_[patchi] - 2*U_.boundaryField()[patchi];
		}else
		{
			pB = p_.boundaryField()[patchi];
			pBB = p_.boundaryField()[patchi];
			rhoB = rho_.boundaryField()[patchi];
			rhoBB = rho_.boundaryField()[patchi];
			UB = 2*U_.boundaryField()[patchi] - U_.boundaryField()[patchi].patchInternalField();
			UBB = 3*UB_[patchi] - 2*U_.boundaryField()[patchi];
		}
	}
}

template<class Type>
void
wasFluxScheme<Type>::sampling
(
	double (exactRiemannSolver::*varFunc)(double)const,
	const exactRiemannSolver& riemann,
	double* var
)const
{
	var[0] = (riemann.*varFunc)(riemann.speed(0)-0.1);
	if(riemann.speed(0)<0.0 && riemann.speed(1)>0.0)
	{
		var[1] = (riemann.*varFunc)(0.0);
	}else
	{
		var[1] = (riemann.*varFunc)(0.5*(riemann.speed(1) + riemann.speed(2)));
	}
	if(riemann.speed(3)<0.0 && riemann.speed(4)>0.0)
	{
		var[2] = (riemann.*varFunc)(0.0);
	}else
	{
		var[2] = (riemann.*varFunc)(0.5*(riemann.speed(2) + riemann.speed(3)));
	}
	var[3] = (riemann.*varFunc)(riemann.speed(4)+0.1);
}

template<class Type>
double
wasFluxScheme<Type>::average
(
	const double* var,
	const double* w
)const
{
	double varAve;

	varAve = 0.5*(var[0] + var[3]);

	for (int i=0; i<3; i++)
	{
		varAve -= w[i]*(var[i+1] - var[i]);
	}

	return varAve;
}

template<class Type>
vector
wasFluxScheme<Type>::averageV
(
	const scalar& speed,
	const vector& vP,
	const vector& vN,
	const tensor& gradVP,
	const tensor& gradVN,
	const vector& d
)const
{
	scalar cfl;
	NVDVTVDV tvd;
	const double dT = p_.mesh().time().deltaTValue();

	cfl = speed*dT/mag(d);
	double r = tvd.r
			(
				cfl,
				vP,
				vN,
				gradVP,
				gradVN,
				d
			);
	double w;
	w = 0.5*sign(cfl)*(1 - (1 - mag(cfl))*limiter(r));
	return 0.5*(vP + vN) - w*(vN - vP);
}

template<class Type>
void
wasFluxScheme<Type>::weight
(
	double* w,
	const double* rho,
	const double* speed,
	const double* rhoL,
	const double* rhoR,
	const vector& d
)const
{
	double cfl;
	const double dT = p_.mesh().time().deltaTValue();

	for (int i=0; i<3; i++)
	{
		cfl = speed[i]*dT/mag(d);
		//TODO to find a better way to define the r function
		double r = rFunc
		        (
		            cfl,
		            rho[i+1] - rho[i],
		            rhoL[i+1] - rhoL[i],
		            rhoR[i+1] - rhoR[i]
		        );
		if(cfl==0.0)
		{
			w[i] = 0.0;
		}else
		{
			w[i] = 0.5*(cfl/::fabs(cfl))*(1 - (1 - ::fabs(cfl))*limiter(r));
		}
	}
}

template<class Type>
double
wasFluxScheme<Type>::limiter
(
	const double r
)const
{
	//SuperA limiter function
	double B;
	if(r<0)
	{
		B = 0.0;
	}else if(r<0.5)
	{
		B = 2*r;
	}else if(r<1)
	{
		B = 1.0;
	}else if(r<2)
	{
		B = r;
	}else if(r<1000)
	{
		B = 2;
	}else
	{
		B = 1;
	}
	return B;
}

template<class Type>
scalar
wasFluxScheme<Type>::rFunc
(
	const scalar& faceFlux,
	const scalar& dPhi,
	const scalar& gradcP,
	const scalar& gradcN
) const
{
	//scalar toll = 1e-6;
	//scalar gradf = phiN - phiP;

	scalar gradcf;

	if (faceFlux > 0)
	{
		gradcf = gradcP;
	}
	else
	{
		gradcf = gradcN;
	}

	if (mag(gradcf) >= 1000*mag(dPhi))
	{
		return 1000*sign(gradcf)*sign(dPhi);
	}
	else
	{
		return gradcf/dPhi;
	}
}

template<class Type>
void
wasFluxScheme<Type>::test
(
)
{
	Info << "Running" << endl;
	riemann();
	Info << "OK" << endl;
}

template<class Type>
wasFluxScheme<Type>::~wasFluxScheme()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
