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

#include "kurganovConvectionScheme.H"
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
kurganovConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& posFaceFlux,
    const surfaceScalarField& negFaceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	surfaceScalarField a_pos(ap_/(ap_ - am_));  //equation.9.b (i.e. Kurganov)
	surfaceScalarField a_neg(1.0 - a_pos);
	surfaceScalarField aSf(am_*a_pos);  //(-1)*equation.10.b (i.e. Kurganov)
	surfaceScalarField aphiv_pos(a_pos*posFaceFlux - aSf);
	surfaceScalarField aphiv_neg(a_neg*negFaceFlux + aSf);

    tmp<GeometricField<Type, fvPatchField, volMesh> > tConvection
    (
        fvc::surfaceIntegrate(interpolate(aphiv_pos, aphiv_neg, vf))
    );

    tConvection().rename
    (
        "convection(" + posFaceFlux.name() + ',' + negFaceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
