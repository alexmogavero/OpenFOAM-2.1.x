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
    Foam::fv::centralConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    centralConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef centralConvectionScheme_H
#define centralConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class centralConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class centralConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

		const surfaceScalarField pos_;

		const surfaceScalarField neg_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        centralConvectionScheme(const centralConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const centralConvectionScheme&);

protected:
        tmp<surfaceInterpolationScheme<Type> > posTinterpScheme_;

        tmp<surfaceInterpolationScheme<Type> > negTinterpScheme_;

        const surfaceScalarField& negFaceFlux_;

public:

    //- Runtime type information
    TypeName("Central");

    // Constructors

        //- Construct from flux and interpolation scheme
        centralConvectionScheme
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
        centralConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& flux,
            Istream& is
        )
        :
			convectionScheme<Type>(mesh, flux),
			pos_(IOobject("pos","dummy",mesh),mesh,dimensionedScalar("pos", dimless, 1.0)),
			neg_(IOobject("neg","dummy",mesh),mesh,dimensionedScalar("neg", dimless, -1.0)),
			posTinterpScheme_(surfaceInterpolationScheme<Type>::New(mesh, pos_, is)),
			negFaceFlux_(static_cast<const objectRegistry&>(mesh).lookupObject<surfaceScalarField>(word(is)))
		{
        	is.rewind();
        	word dummy(is);
        	negTinterpScheme_ = surfaceInterpolationScheme<Type>::New(mesh, neg_, is);
        	word dummy1(is);
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

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
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
#   include "centralConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
