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

#include "goudonovFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type>
void
goudonovFluxScheme<Type>::calculate
(
	surfaceScalarField& pAve,
	surfaceScalarField& rhoAve,
	surfaceVectorField& UAve
)
{
	const fvMesh& mesh = p_.mesh();
	const labelUList& owner = mesh.owner();
	const labelUList& neighbour = mesh.neighbour();

	forAll(pAve,face)
	{
		label own = owner[face];
		label nei = neighbour[face];

		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;
		const vector vL = U_[own] - uL*n;
		const vector vR = U_[nei] - uR*n;

		//TODO implement gamma from thermal model
		/*if(face==308 || face==68)
		{
			exactRiemannSolver::debug = true;
		}else
		{
			exactRiemannSolver::debug = false;
		}*/
		exactRiemannSolver riemann
		        (
		        	p_[own],
		        	p_[nei],
		        	rho_[own],
		        	rho_[nei],
		        	uL,
		        	uR,
		        	1.4
		        );
		riemann.solve();

		pAve[face] = riemann.p(0);
		rhoAve[face] = riemann.rho(0);
		double uAve = riemann.u(0);
		vector vAve = vR + pos(riemann.speed(2))*(vL - vR);  //first order upwind with the star velocity
		UAve[face] = uAve*n + vAve;

		/*if(face==308 || face==68)
		{
			Info << face << " " << own << "-" << nei << endl;
			Info << pAve[face] << " " << rhoAve[face] << " " << uAve << endl;
			Info << p_[own] << " " << p_[nei] << endl;
		}*/
	}

	surfaceScalarField::GeometricBoundaryField& bP = pAve.boundaryField();
	forAll(bP, patchi)
	{
		const scalarField pL = p_.boundaryField()[patchi].patchInternalField();
		const scalarField pR = p_.boundaryField()[patchi];
		const scalarField rhoL = rho_.boundaryField()[patchi].patchInternalField();
		const scalarField rhoR = rho_.boundaryField()[patchi];
		const vectorField UL = U_.boundaryField()[patchi].patchInternalField();
		const vectorField UR = 2*U_.boundaryField()[patchi] - UL;
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

			//TODO implement gamma from thermal model
			exactRiemannSolver riemann(pL[face],pR[face],rhoL[face],rhoR[face],uL,uR,1.4);
			riemann.solve();

			pAveB[face] = riemann.p(0);
			rhoAveB[face] = riemann.rho(0);
			double uAve = riemann.u(0);

			vector vAve = vR + pos(riemann.speed(2))*(vL - vR);
			UAveB[face] = uAve*n + vAve;
		}
	}
}

template<class Type>
goudonovFluxScheme<Type>::~goudonovFluxScheme()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
