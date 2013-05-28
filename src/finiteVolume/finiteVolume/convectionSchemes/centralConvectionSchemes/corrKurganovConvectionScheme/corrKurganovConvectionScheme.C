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

#include "corrKurganovConvectionScheme.H"
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
tmp<GeometricField<Type, fvPatchField, volMesh> >
corrKurganovConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > uPos(this->posTinterpScheme_().interpolate(vf));
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > uNeg(this->negTinterpScheme_().interpolate(vf));
	surfaceScalarField a_pos(this->ap_/(this->ap_ - this->am_));  //equation.9.b (i.e. Kurganov)
	surfaceScalarField a_neg(1.0 - a_pos);
	surfaceScalarField aSf(this->am_*a_pos);  //(-1)*equation.10.b (i.e. Kurganov)
	surfaceScalarField aphiv_pos(a_pos*posFaceFlux - aSf);
	surfaceScalarField aphiv_neg(a_neg*negFaceFlux + aSf);

    GeometricField<Type, fvsPatchField, surfaceMesh> flux
    (
    		aphiv_pos*uPos() + aphiv_neg*uNeg() - (this->ap_*this->am_)*correction(posFaceFlux,negFaceFlux,uPos(),uNeg())
    );

    tmp<GeometricField<Type, fvPatchField, volMesh> > tConvection
	(
			fvc::surfaceIntegrate(flux)
	);

    tConvection().rename
    (
        "convection(" + posFaceFlux.name() + ',' + negFaceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
corrKurganovConvectionScheme<Type>::wIntermediate
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& uPos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& uNeg
) const
{
	return	((this->ap_ - negFaceFlux)*uNeg - (this->am_ - posFaceFlux)*uPos)/(this->ap_ - this->am_);
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
corrKurganovConvectionScheme<Type>::correction
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& uPos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& uNeg
) const
{
	GeometricField<Type, fvsPatchField, surfaceMesh> wInt(wIntermediate(posFaceFlux,negFaceFlux,uPos,uNeg));
	GeometricField<Type, fvsPatchField, surfaceMesh> diffPos((wInt - uPos)/(this->ap_ - this->am_));
	GeometricField<Type, fvsPatchField, surfaceMesh> diffNeg((uNeg - wInt)/(this->ap_ - this->am_));

	GeometricField<Type, fvsPatchField, surfaceMesh> minDiff(min(diffPos,diffNeg));
	return	max(minDiff,0.0*minDiff) + min(minDiff,0.0*minDiff);  //minMod function
}

/*template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
corrKurganovConvectionScheme<Type>::minMod
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& c1,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& c2
) const
{
	dimensioned<Type> zero("zero", c1.dimensions(), Type::zero);

	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > cPos(min(c1,c2));
	c

	return min(c1,c2)().max(zero) + min(c1,c2)().min(zero);
}

template<>
tmp<GeometricField<double, fvsPatchField, surfaceMesh> >
corrKurganovConvectionScheme<double>::minMod
(
    const GeometricField<double, fvsPatchField, surfaceMesh>& diffPos,
    const GeometricField<double, fvsPatchField, surfaceMesh>& diffNeg
) const
{
	dimensioned<double> zero("zero", diffPos.dimensions(), 0.0);

	return min(diffPos,diffNeg)().max(zero) + min(diffPos,diffNeg)().min(zero);
}*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
