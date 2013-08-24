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

#include "corrKurganovFluxScheme.H"
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

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::fluxT
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > uPos(this->interpolatePos(vf));
	tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > uNeg(this->interpolateNeg(vf));

	return this->aphiv_pos_*uPos() + this->aphiv_neg_*uNeg() -
			(this->ap_*this->am_)*correction(this->phiv_pos_,this->phiv_pos_,uPos(),uNeg());
}

template<class Type>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::flux
(
	const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
	return fluxT(vf);
}

template<class Type>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::flux
(
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	return fluxT(vf);
}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> >
corrKurganovFluxScheme<Type>::wIntermediate
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type1, fvsPatchField, surfaceMesh>& uPos,
    const GeometricField<Type1, fvsPatchField, surfaceMesh>& uNeg
) const
{
	return	((this->ap_ - negFaceFlux)*uNeg - (this->am_ - posFaceFlux)*uPos)/(this->ap_ - this->am_);
}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> >
corrKurganovFluxScheme<Type>::correction
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type1, fvsPatchField, surfaceMesh>& uPos,
    const GeometricField<Type1, fvsPatchField, surfaceMesh>& uNeg
) const
{
	GeometricField<Type1, fvsPatchField, surfaceMesh> wInt(wIntermediate(posFaceFlux,negFaceFlux,uPos,uNeg));
	GeometricField<Type1, fvsPatchField, surfaceMesh> diffPos((wInt - uPos)/(this->ap_ - this->am_));
	GeometricField<Type1, fvsPatchField, surfaceMesh> diffNeg((uNeg - wInt)/(this->ap_ - this->am_));

	GeometricField<Type1, fvsPatchField, surfaceMesh> minDiff(min(diffPos,diffNeg));
	return	max(minDiff,0.0*minDiff) + min(minDiff,0.0*minDiff);  //minMod function
}

template<class Type>
template<class Type1>
tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::diffusionT
(
	const GeometricField<Type1, fvPatchField, volMesh>& vf
) const
{
	const tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > uPos(this->interpolatePos(vf));
	const tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > uNeg(this->interpolateNeg(vf));

	return this->aSf_*(uNeg() - uPos()) -
			(this->ap_*this->am_)*correction(this->phiv_pos_,this->phiv_pos_,uPos(),uNeg());
}

template<class Type>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::diffusion
(
	const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
	return diffusionT(vf);
}

template<class Type>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > corrKurganovFluxScheme<Type>::diffusion
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
