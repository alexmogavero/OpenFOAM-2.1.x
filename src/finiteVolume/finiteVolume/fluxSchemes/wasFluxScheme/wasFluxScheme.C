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
wasFluxScheme<Type>::calculateOld
(
	surfaceScalarField& pAve,
	surfaceScalarField& rhoAve,
	surfaceVectorField& UAve
)
{
	const fvMesh& mesh = p_.mesh();
	const labelUList& owner = mesh.owner();
	const labelUList& neighbour = mesh.neighbour();
	const volVectorField& C = mesh.C();
	const volVectorField gradRho(fvc::grad(rho_));
	const volVectorField gradP(fvc::grad(p_));
	const volTensorField gradU(fvc::grad(U_));

	forAll(pAve,face)
	{
		label own = owner[face];
		label nei = neighbour[face];

		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;
		const vector vL = U_[own] - uL*n;
		const vector vR = U_[nei] - uR*n;
		const vector gradUL = gradU[own] & n;
		const vector gradUR = gradU[nei] & n;
		const tensor gradVL = gradU[own] - (gradU[own] & n)*n;
		const tensor gradVR = gradU[nei] - (gradU[nei] & n)*n;

		//TODO implement gamma from thermal model
		exactRiemannSolver riemann(p_[own],p_[nei],rho_[own],rho_[nei],uL,uR,1.4);
		riemann.solve();

		const vector d = C[nei] - C[own];
		Info << own << "-" << nei << endl;
		Info << p_[nei] << " " << gradP[own] << endl;
		Info << p_[nei]-2*(gradP[own] & d) << " " << p_[own] << " " << rho_[nei]-2*(gradRho[own] & d) << " " << rho_[own] << " " <<
				uR-2*(gradUL & d) << " " << uL << endl;
		exactRiemannSolver riemannL(p_[nei]-2*(gradP[own] & d),p_[own],
				rho_[nei]-2*(gradRho[own] & d),rho_[own],
				uR-2*(gradUL & d),uL,1.4);
		riemannL.solve();
		Info << p_[nei] << " " << p_[own]+2*(gradP[nei] & d) << " " << rho_[nei] << " " << rho_[own]+2*(gradRho[nei] & d) << " " <<
				uR << " " << uL+2*(gradUR & d) << endl;
		exactRiemannSolver riemannR(p_[nei],p_[own]+2*(gradP[nei] & d),
				rho_[nei],rho_[own]+2*(gradRho[nei] & d),
				uR,uL+2*(gradUR & d),1.4);
		riemannR.solve();

		double w[3];
		//weight(w,riemann,gradRho[own],gradRho[nei],d);
		weight1(w,riemann,riemannL,riemannR,d);

		pAve[face] = average(&exactRiemannSolver::p,riemann,w);
		rhoAve[face] = average(&exactRiemannSolver::rho,riemann,w);
		double uAve = average(&exactRiemannSolver::u,riemann,w);
		vector vAve = averageV(riemann,vL,vR,gradVL,gradVR,d);
		UAve[face] = uAve*n + vAve;
	}

	const volScalarField::GeometricBoundaryField& bP = p_.boundaryField();
	forAll(bP, patchi)
	{
		const scalarField pL = bP[patchi].patchInternalField();
		const scalarField pR = bP[patchi];
		const scalarField rhoL = rho_.boundaryField()[patchi].patchInternalField();
		const scalarField rhoR = rho_.boundaryField()[patchi];
		const vectorField UL = U_.boundaryField()[patchi].patchInternalField();
		const vectorField UR = 2*U_.boundaryField()[patchi] - UL;
		const tensorField gradUL = gradU.boundaryField()[patchi].patchInternalField();
		const tensorField gradUR = 2*gradU.boundaryField()[patchi] - gradU.boundaryField()[patchi].patchInternalField();
		const vectorField gradRhoL = gradRho.boundaryField()[patchi].patchInternalField();
		const vectorField gradRhoR = 2*gradRho.boundaryField()[patchi] - gradRhoL;
		const vectorField gradPL = gradP.boundaryField()[patchi].patchInternalField();
		const vectorField gradPR = 2*gradP.boundaryField()[patchi] - gradPL;
		const vectorField d = 2*mesh.boundary()[patchi].delta();
		scalarField& pAveB = pAve.boundaryField()[patchi];
		scalarField& rhoAveB = rhoAve.boundaryField()[patchi];
		vectorField& UAveB = UAve.boundaryField()[patchi];

		forAll(pAveB,face)
		{
			const vector n = mesh.boundary()[patchi].Sf()[face]/mesh.boundary()[patchi].magSf()[face];
			const double uL = UL[face] & n;
			const double uR = UR[face] & n;
			const vector vL = UL[face] - uL*n;
			const vector vR = UR[face] - uR*n;
			const vector graduL = gradUL[face] & n;
			const vector graduR = gradUR[face] & n;
			const tensor gradVL = gradUL[face] - (gradUL[face] & n)*n;
			const tensor gradVR = gradUR[face] - (gradUR[face] & n)*n;

			//TODO implement gamma from thermal model
			Info << mesh.boundary()[patchi].name() << endl;
			Info << pL[face] << " " << pR[face] << " " << rhoL[face] << " " << rhoR[face] << " " <<
					uL << " " << uR << endl;
			exactRiemannSolver riemann(pL[face],pR[face],rhoL[face],rhoR[face],uL,uR,1.4);
			riemann.solve();

			scalar pLL(pR[face]-2*(gradPL[face] & d[face]));
			scalar rhoLL(rhoR[face]-2*(gradRhoL[face] & d[face]));
			scalar uLL(uR-2*(graduL & d[face]));
			Info << gradPL[face] << " " << d[face] << endl;
			Info << pLL << " " << pL[face] << " " << rhoLL << " " << rhoL[face] << " " <<
					uLL << " " << uL << endl;
			exactRiemannSolver riemannL(pLL,pL[face],rhoLL,rhoL[face],uLL,uL,1.4);
			riemannL.solve();
			scalar pRR(pR[face]); //pRR(pL[face]+2*(gradPR[face] & d[face]));
			scalar rhoRR(rhoR[face]); //(rhoL[face]+2*(gradRhoR[face] & d[face]));
			scalar uRR(uL+2*(graduR & d[face]));
			Info << pR[face] << " " << pRR << " " << rhoR[face] << " " << rhoRR << " " <<
					uR << " " << uRR << endl;
			exactRiemannSolver riemannR(pR[face],pRR,rhoR[face],rhoRR,uR,uRR,1.4);
			riemannR.solve();

			double w[3];
			//weight(w,riemann,gradRhoL[face],gradRhoR[face],2*d[face]);
			weight1(w,riemann,riemannL,riemannR,d[face]);

			pAveB[face] = average(&exactRiemannSolver::p,riemann,w);
			rhoAveB[face] = average(&exactRiemannSolver::rho,riemann,w);
			double uAve = average(&exactRiemannSolver::u,riemann,w);

			vector vAve = averageV(riemann,vL,vR,gradVL,gradVR,d[face]);
			UAveB[face] = uAve*n + vAve;
		}
	}
}

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
	const volVectorField gradRho(fvc::grad(rho_));
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
		//const double rho[4] = {rho_[own], rho1_[face], rho2_[face], rho_[nei]};
		const double rho[4] = {rho0_[face], rho1_[face], rho2_[face], rho3_[face]};
		const double speed[3] = {speed0_[face], speed1_[face], speed2_[face]};

		/*const double rhoL[4] = {rho_[nei]-2*(d & gradRho[own]), rho1_[face]-(d & gradRho1[own]),
				rho2_[face]-(d & gradRho2[own]), rho_[own]};
		const double rhoR[4] = {rho_[nei], rho1_[face]+(d & gradRho1[nei]),
				rho2_[face]+(d & gradRho2[nei]), rho_[own]+2*(d & gradRho[nei])};*/
		const double rhoL[4] = {rho0_[face]-(d & gradRho0[own]), rho1_[face]-(d & gradRho1[own]),
				rho2_[face]-(d & gradRho2[own]), rho3_[face]-(d & gradRho3[own])};
		const double rhoR[4] = {rho0_[face]+(d & gradRho0[nei]), rho1_[face]+(d & gradRho1[nei]),
				rho2_[face]+(d & gradRho2[nei]), rho3_[face]+(d & gradRho3[nei])};

		/*Info << face << " " << own << " " << nei << endl;
		Info << "rho  " << rho[0] << " " << rho[1] << " " << rho[2] << " " << rho[3] << endl;
		Info << "rhoL " << rhoL[0] << " " << rhoL[1] << " " << rhoL[2] << " " << rhoL[3] << endl;
		Info << "rhoR " << rhoR[0] << " " << rhoR[1] << " " << rhoR[2] << " " << rhoR[3] << endl;*/
		weight(w, rho, speed, rhoL, rhoR, d);

		//const double p[4] = {p_[own], p1_[face], p2_[face], p_[nei]};
		const double p[4] = {p0_[face], p1_[face], p2_[face], p3_[face]};
		pAve[face] = average(p,w);
		rhoAve[face] = average(rho,w);
		//const double u[4] = {uL, u1_[face], u2_[face], uR};
		const double u[4] = {u0_[face], u1_[face], u2_[face], u3_[face]};
		const scalar uAve = average(u,w);

		vector vAve = averageV(speed1_[face],vL,vR,gradVL,gradVR,d);
		UAve[face] = uAve*n + vAve;

		/*if(face==237)
		{
			Info << face << endl;
			Info << speed[0] << " " << speed[1] << " " << speed[2] << endl;
			Info << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << " " << uAve << endl;
			Info << w[0] << " " << w[1] << " " << w[2] << endl;
			Info << rho[0] << " " << rho[1] << " " << rho[2] << " " << rho[3] << endl;
			Info << rhoL[0] << " " << rhoL[1] << " " << rhoL[2] << " " << rhoL[3] << endl;
			Info << rhoR[0] << " " << rhoR[1] << " " << rhoR[2] << " " << rhoR[3] << endl;
		}*/
	}

	const surfaceScalarField::GeometricBoundaryField& bP = p1_.boundaryField();
	forAll(bP, patchi)
	{
		const scalarField pL = p_.boundaryField()[patchi].patchInternalField();
		const scalarField pR = p_.boundaryField()[patchi];
		const scalarField p0 = p0_.boundaryField()[patchi];
		const scalarField p1 = p1_.boundaryField()[patchi];
		const scalarField p2 = p2_.boundaryField()[patchi];
		const scalarField p3 = p3_.boundaryField()[patchi];
		const scalarField rhoL = rho_.boundaryField()[patchi].patchInternalField();
		const scalarField rhoR = rho_.boundaryField()[patchi];
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
		const vectorField gradRhoL = gradRho.boundaryField()[patchi].patchInternalField();
		const vectorField gradRhoR = 2*gradRho.boundaryField()[patchi] - gradRhoL;
		const vectorField gradRho0L = gradRho0.boundaryField()[patchi].patchInternalField();
		const vectorField gradRho1L = gradRho1.boundaryField()[patchi].patchInternalField();
		const vectorField gradRho1R = 2*gradRho1.boundaryField()[patchi] - gradRho1L;
		const vectorField gradRho2L = gradRho2.boundaryField()[patchi].patchInternalField();
		const vectorField gradRho2R = 2*gradRho2.boundaryField()[patchi] - gradRho2L;
		const vectorField gradRho3L = gradRho3.boundaryField()[patchi].patchInternalField();

		const vectorField d = 2*mesh.boundary()[patchi].delta();

		scalarField& pAveB = pAve.boundaryField()[patchi];
		scalarField& rhoAveB = rhoAve.boundaryField()[patchi];
		vectorField& UAveB = UAve.boundaryField()[patchi];

		forAll(pAveB,face)
		{
			double w[3];
			//const double rhoB[4] = {rhoL[face], rho1[face], rho2[face], rhoR[face]};
			const double rhoB[4] = {rho0[face], rho1[face], rho2[face], rho3[face]};
			const double speedB[3] = {speed0[face], speed1[face], speed2[face]};

			/*const double rhoLB[4] = {rhoR[face]-2*(d[face] & gradRhoL[face]), rho1[face]-(d[face] & gradRho1L[face]),
					rho2[face]-(d[face] & gradRho2L[face]), rhoL[face]};
			const double rhoRB[4] = {rhoB[0],rhoB[1],rhoB[2],rhoB[3]}; //{rhoR[face], rho1[face]+(d[face] & gradRho1R[face]),
					//rho2[face]+(d[face] & gradRho2R[face]), rhoL[face]+2*(d[face] & gradRhoR[face])};*/
			const double rhoLB[4] = {rho0[face]-(d[face] & gradRho0L[face]), rho1[face]-(d[face] & gradRho1L[face]),
					rho2[face]-(d[face] & gradRho2L[face]), rho3[face]-(d[face] & gradRho3L[face])};
			const double rhoRB[4] = {rho0B[face],rho1B[face],rho2B[face],rho3B[face]};

			weight(w, rhoB, speedB, rhoLB, rhoRB, d[face]);

			//const double pB[4] = {pL[face], p1[face], p2[face], pR[face]};
			const double pB[4] = {p0[face], p1[face], p2[face], p3[face]};
			pAveB[face] = average(pB,w);
			rhoAveB[face] = average(rhoB,w);
			const vector n = mesh.boundary()[patchi].Sf()[face]/mesh.boundary()[patchi].magSf()[face];
			const double uL = UL[face] & n;
			const double uR = UR[face] & n;
			const vector vL = UL[face] - uL*n;
			const vector vR = UR[face] - uR*n;
			const tensor gradVL = gradUL[face] - (gradUR[face] & n)*n;
			const tensor gradVR = gradUR[face] - (gradUR[face] & n)*n;
			//const double uB[4] = {uL, u1[face], u2[face], uR};
			const double uB[4] = {u0[face], u1[face], u2[face], u3[face]};
			const scalar uAve = average(uB,w);

			vector vAve = averageV(speed1[face],vL,vR,gradVL,gradVR,d[face]);
			UAveB[face] = uAve*n + vAve;

			/*if((mesh.boundary()[patchi].name()=="wall") && face==0)
			{
				Info << face << " wall" << endl;
				Info << speedB[0] << " " << speedB[1] << " " << speedB[2] << endl;
				Info << uB[0] << " " << uB[1] << " " << uB[2] << " " << uB[3] << " " << uAve << endl;
				Info << pB[0] << " " << pB[1] << " " << pB[2] << " " << pB[3] << " " << pAveB[face] << endl;
				Info << rhoB[0] << " " << rhoB[1] << " " << rhoB[2] << " " << rhoB[3] << " " << pAveB[face] << endl;
				Info << w[0] << " " << w[1] << " " << w[2] << endl;
				Info << rhoB[0] << " " << rhoB[1] << " " << rhoB[2] << " " << rhoB[3] << endl;
				Info << rhoLB[0] << " " << rhoLB[1] << " " << rhoLB[2] << " " << rhoLB[3] << endl;
				Info << rhoRB[0] << " " << rhoRB[1] << " " << rhoRB[2] << " " << rhoRB[3] << endl;
			}*/
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

	forAll(p1_,face)
	{
		label own = owner[face];
		label nei = neighbour[face];

		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;

		//TODO implement gamma from thermal model
		exactRiemannSolver riemann(p_[own],p_[nei],rho_[own],rho_[nei],uL,uR,1.4);
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
		const scalarField pR = pB_[patchi]; //p_.boundaryField()[patchi];
		const scalarField pRR = pBB_[patchi]; //p_.boundaryField()[patchi];
		const scalarField rhoL = rho_.boundaryField()[patchi].patchInternalField();
		const scalarField rhoR = rhoB_[patchi]; //rho_.boundaryField()[patchi];
		const scalarField rhoRR = rhoBB_[patchi]; //rho_.boundaryField()[patchi];
		const vectorField UL = U_.boundaryField()[patchi].patchInternalField();
		const vectorField UR = UB_[patchi]; //2*U_.boundaryField()[patchi] - UL;
		const vectorField URR = UBB_[patchi]; //3*UR - 2*U_.boundaryField()[patchi];
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

			//TODO implement gamma from thermal model
			exactRiemannSolver riemann(pL[face],pR[face],rhoL[face],rhoR[face],uL,uR,1.4);
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

			//TODO implement gamma from thermal model
			exactRiemannSolver riemannR(pR[face],pRR[face],rhoR[face],rhoRR[face],uR,uRR,1.4);
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
	double (exactRiemannSolver::*varFunc)(double)const,
	const exactRiemannSolver& riemann,
	const double* w
)const
{
	double var[4];
	double varAve;

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

	varAve = 0.5*(var[0] + var[3]);

	for (int i=0; i<3; i++)
	{
		varAve -= w[i]*(var[i+1] - var[i]);
	}

	return varAve;
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
	const exactRiemannSolver& riemann,
	const vector& vP,
	const vector& vN,
	const tensor& gradVP,
	const tensor& gradVN,
	const vector& d
)const
{
	double speed;
	double cfl;
	NVDVTVDV tvd;
	const double dT = p_.mesh().time().deltaTValue();

	speed = riemann.speed(2);

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
	if(cfl==0.0)
	{
		w = 0.0;
	}else
	{
		w = 0.5*(cfl/::fabs(cfl))*(1 - (1 - ::fabs(cfl))*limiter(r));
	}
	return 0.5*(vP + vN) - w*(vN - vP);
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
	/*if(cfl==0.0)
	{
		w = 0.0;
	}else
	{
		w = 0.5*(cfl/::fabs(cfl))*(1 - (1 - ::fabs(cfl))*limiter(r));
	}*/
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
void
wasFluxScheme<Type>::weight1
(
	double* w,
	const exactRiemannSolver& riemann,
	const exactRiemannSolver& riemannL,
	const exactRiemannSolver& riemannR,
	const vector& d
)const
{
	scalar speed[3];
	scalar cfl;
	//NVDTVD tvd;
	const scalar dT = p_.mesh().time().deltaTValue();
	scalar rho[4];
	scalar rhoL[4];
	scalar rhoR[4];
	const scalar toll=1e-6;

	speed[0] = riemann.speed(0);
	speed[1] = riemann.speed(2);
	speed[2] = riemann.speed(4);
	if(riemann.speed(0)<0.0 && riemann.speed(1)>0.0)
	{
		//Info << "sonic" << endl;
	}
	if(riemann.speed(3)*riemann.speed(4)<0.0)
	{
		//Info << "sonic" << endl;
	}
	rho[0] = riemann.rho(riemann.speed(0)-0.1);
	if(riemann.speed(0)<0.0 && riemann.speed(1)>0.0)
	{
		rho[1] = riemann.rho(0.0);
	}else
	{
		rho[1] = riemann.rho(0.5*(riemann.speed(1) + riemann.speed(2)));
	}
	if(riemann.speed(3)<0.0 && riemann.speed(4)>0.0)
	{
		rho[2] = riemann.rho(0.0);
	}else
	{
		rho[2] = riemann.rho(0.5*(riemann.speed(2) + riemann.speed(3)));
	}
	rho[3] = riemann.rho(riemann.speed(3)+0.1);
	rhoL[0] = riemannL.rho(riemannL.speed(0)-0.1);
	if(riemannL.speed(0)<0.0 && riemannL.speed(1)>0.0)
	{
		rhoL[1] = riemannL.rho(0.0);
	}else
	{
		rhoL[1] = riemannL.rho(0.5*(riemannL.speed(1) + riemannL.speed(2)));
	}
	if(riemannL.speed(3)<0.0 && riemannL.speed(4)>0.0)
	{
		rhoL[2] = riemannL.rho(0.0);
	}else
	{
		rhoL[2] = riemannL.rho(0.5*(riemannL.speed(2) + riemannL.speed(3)));
	}
	rhoL[3] = riemannL.rho(riemannL.speed(3)+0.1);
	rhoR[0] = riemannR.rho(riemannR.speed(0)-0.1);
	if(riemannR.speed(0)<0.0 && riemannR.speed(1)>0.0)
	{
		rhoR[1] = riemannR.rho(0.0);
	}else
	{
		rhoR[1] = riemannR.rho(0.5*(riemannR.speed(1) + riemannR.speed(2)));
	}
	if(riemannR.speed(3)<0.0 && riemannR.speed(4)>0.0)
	{
		rhoR[2] = riemannR.rho(0.0);
	}else
	{
		rhoR[2] = riemannR.rho(0.5*(riemannR.speed(2) + riemannR.speed(3)));
	}
	rhoR[3] = riemannR.rho(riemannR.speed(3)+0.1);

	for (int i=0; i<3; i++)
	{
		cfl = speed[i]*dT/mag(d);
		scalar r;
		if(cfl>=0)
		{
			scalar diffLocal = sign(rho[i+1] - rho[i])*max(mag(rho[i+1] - rho[i]), toll*0.5*(rho[i+1] + rho[i]));
			scalar diffUp = sign(rhoL[i+1] - rhoL[i])*max(mag(rhoL[i+1] - rhoL[i]), toll*0.5*(rhoL[i+1] + rhoL[i]));
			r = diffUp/diffLocal;
		}else
		{
			scalar diffLocal = sign(rho[i+1] - rho[i])*max(mag(rho[i+1] - rho[i]), toll*0.5*(rho[i+1] + rho[i]));
			scalar diffUp = sign(rhoR[i+1] - rhoR[i])*max(mag(rhoR[i+1] - rhoR[i]), toll*0.5*(rhoR[i+1] + rhoR[i]));
			r = diffUp/diffLocal;
		}
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
