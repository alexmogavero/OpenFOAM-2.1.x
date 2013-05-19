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

Class
    Foam::fv::kurganovConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    kurganovConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef kurganovConvectionScheme_H
#define kurganovConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class kurganovConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class kurganovConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

		Istream& is_;

		const tmp<surfaceScalarField> pos_;

		const tmp<surfaceScalarField> neg_;

        tmp<surfaceInterpolationScheme<Type> > posTinterpScheme_;

        tmp<surfaceInterpolationScheme<Type> > negTinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kurganovConvectionScheme(const kurganovConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const kurganovConvectionScheme&);


public:

    //- Runtime type information
    TypeName("Kurganov");


    // Constructors

        //- Construct from flux and interpolation scheme
        kurganovConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type> >& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            posTinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        kurganovConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& pos,
            Istream& is
        )
        :
			convectionScheme<Type>(mesh, pos),
			is_(is),
			pos_(pos),
			neg_(-pos),
			posTinterpScheme_(surfaceInterpolationScheme<Type>::New(mesh, pos_, is)),
			negTinterpScheme_(surfaceInterpolationScheme<Type>::New(mesh, neg_, is.rewind()))
		{
        	posTinterpScheme_().extrapolate_ = true;
        	negTinterpScheme_().extrapolate_ = true;
		}


    // Member Functions

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
        (
            const surfaceScalarField&,
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kurganovConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
