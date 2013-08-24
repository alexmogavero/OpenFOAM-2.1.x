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

#include "kurganovFluxScheme.H"
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

/*template<class Type>
void kurganovFluxScheme<Type>::phiv
(
)
{
	const fvMesh& mesh(this->mesh());


}*/

template<class Type>
void kurganovFluxScheme<Type>::weight
(
)
{
	const fvMesh& mesh(this->mesh());
	const surfaceVectorField U_pos(this->interpolatePos(this->U_));
	const surfaceVectorField U_neg(this->interpolateNeg(this->U_));
	phiv_pos_ = (U_pos & mesh.Sf());
	phiv_neg_ = (U_neg & mesh.Sf());

	const volScalarField& psi = this->thermo_.psi();
	volScalarField rPsi(1.0/psi);
	volScalarField c(sqrt(this->thermo_.Cp()/this->thermo_.Cv()*rPsi));

	surfaceScalarField cSf_pos(this->interpolatePos(c)*mesh.magSf());
	surfaceScalarField cSf_neg(this->interpolateNeg(c)*mesh.magSf());

	dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
	ap_ = max(max(phiv_pos_ + cSf_pos, phiv_neg_ + cSf_neg), v_zero);  //equation.8.a
	am_ = min(min(phiv_pos_ - cSf_pos, phiv_neg_ - cSf_neg), v_zero);  //(-1)*equation.8.b

	a_pos_ = ap_/(ap_ - am_);  //equation.9.b (i.e. Kurganov)
	a_neg_ = 1.0 - a_pos_;
	aSf_ = am_*a_pos_;  //(-1)*equation.10.b (i.e. Kurganov)
	aphiv_pos_ = a_pos_*phiv_pos_ - aSf_;
	aphiv_neg_ = a_neg_*phiv_neg_ + aSf_;

}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::fluxT
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	return aphiv_pos_*this->interpolatePos(vf) + aphiv_neg_*this->interpolateNeg(vf);
}

template<class Type>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::flux
(
	const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
	return fluxT(vf);
}

template<class Type>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::flux
(
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	return fluxT(vf);
}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::interpolateT
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	return a_pos_*this->interpolatePos(vf) + a_neg_*this->interpolateNeg(vf);
}

template<class Type>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::interpolate
(
	const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
	return interpolateT(vf);
}

template<class Type>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::interpolate
(
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	return interpolateT(vf);
}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::diffusionT
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	return aSf_*(this->interpolateNeg(vf) - this->interpolatePos(vf));
}

template<class Type>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::diffusion
(
	const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
	return diffusionT(vf);
}

template<class Type>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > kurganovFluxScheme<Type>::diffusion
(
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	return diffusionT(vf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
