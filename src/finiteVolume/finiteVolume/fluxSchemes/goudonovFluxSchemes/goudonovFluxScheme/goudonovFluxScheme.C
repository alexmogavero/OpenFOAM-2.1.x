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
	const tmp<volScalarField> Cp = thermo_.Cp();
	const tmp<volScalarField> Cv = thermo_.Cv();
	const surfaceScalarField gamma(fvc::interpolate(Cp()/Cv()));

	forAll(pAve,face)
	{
		label own = owner[face];
		label nei = neighbour[face];

		const vector n = mesh.Sf()[face]/mesh.magSf()[face];
		const double uL = U_[own] & n;
		const double uR = U_[nei] & n;
		const vector vL = U_[own] - uL*n;
		const vector vR = U_[nei] - uR*n;

		exactRiemannSolver riemann
		        (
		        	p_[own],
		        	p_[nei],
		        	rho_[own],
		        	rho_[nei],
		        	uL,
		        	uR,
		        	gamma[face]
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
		const scalarField gammab = gamma.boundaryField()[patchi];
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

			exactRiemannSolver riemann(pL[face],pR[face],rhoL[face],rhoR[face],uL,uR,gammab[face]);
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
void goudonovFluxScheme<Type>::flux
(
	surfaceScalarField& rhoFlux,
	surfaceVectorField& UFlux,
	surfaceScalarField& EFlux
)
{
	const fvMesh& mesh(this->mesh());
	surfaceScalarField pAve
	(
		IOobject
		(
			"pAve",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("pAve", p_.dimensions(), 0.0)
	);
	surfaceScalarField rhoAve
	(
		IOobject
		(
			"rhoAve",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("rhoAve", rho_.dimensions(), 0.0)
	);
	surfaceVectorField UAve
	(
		IOobject
		(
			"UAve",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedVector("UAve", U_.dimensions(), vector::zero)
	);

	calculate(pAve,rhoAve,UAve);

	rhoFlux = rhoAve*(UAve & mesh.Sf());

	UFlux = rhoFlux*UAve + pAve*mesh.Sf();

	surfaceScalarField R(fvc::interpolate(thermo_.Cp() - thermo_.Cv()));
	surfaceScalarField TAve(pAve/(R*rhoAve));
	const labelList& own = mesh.faceOwner();
	surfaceScalarField eAve
	(
		IOobject
		(
			"eAve",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("eAve", R.dimensions()*TAve.dimensions(), 0.0)
	);
	eAve.internalField() = thermo_.e(TAve.internalField(),own);
	forAll(TAve.boundaryField(),patchi)
	{
		eAve.boundaryField()[patchi] = thermo_.e(TAve.boundaryField()[patchi],patchi);
	}
	EFlux = (UAve & mesh.Sf())*(eAve*rhoAve + 0.5*rhoAve*magSqr(UAve) + pAve);
}

template<class Type>
goudonovFluxScheme<Type>::~goudonovFluxScheme()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
