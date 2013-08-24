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

#include "centralFluxScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
template <class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > centralFluxScheme<Type>::interpolatePos
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	is_.rewind();
	word dummy(is_);
	const fvMesh& mesh(this->mesh());
	tmp<surfaceInterpolationScheme<Type1> > interpScheme(surfaceInterpolationScheme<Type1>::New(mesh, pos_, is_));
	interpScheme().extrapolate_ = true;

	return interpScheme().interpolate(vf);
}

template <class Type>
template <class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > centralFluxScheme<Type>::interpolateNeg
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	is_.rewind();
	word dummy(is_);
	const fvMesh& mesh(this->mesh());
	tmp<surfaceInterpolationScheme<Type1> > interpScheme(surfaceInterpolationScheme<Type1>::New(mesh, neg_, is_));
	interpScheme().extrapolate_ = true;

	return interpScheme().interpolate(vf);
}

template<class Type>
void centralFluxScheme<Type>::flux
(
	surfaceScalarField& rhoFlux,
	surfaceVectorField& UFlux,
	surfaceScalarField& EFlux
)
{
	weight();

	const volScalarField& p(thermo_.p());
	const volScalarField& e(thermo_.e());
	const fvMesh& mesh(this->mesh());

	rhoFlux = flux(rho_);

	volVectorField rhoU(rho_*U_);
	UFlux = flux(rhoU) + interpolate(p)*mesh.Sf();

	volScalarField E(rho_*(e + 0.5*magSqr(U_)) + p);
	EFlux = flux(E) - diffusion(p);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
