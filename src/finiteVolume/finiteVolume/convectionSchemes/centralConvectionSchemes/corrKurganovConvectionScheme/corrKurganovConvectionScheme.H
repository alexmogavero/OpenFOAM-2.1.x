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
    Foam::fv::corrKurganovConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    corrKurganovConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef corrKurganovConvectionScheme_H
#define corrKurganovConvectionScheme_H

#include "kurganovConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class corrKurganovConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class corrKurganovConvectionScheme
:
    public fv::kurganovConvectionScheme<Type>
{
    // Private data

        //const surfaceScalarField& ap_;

        //const surfaceScalarField& am_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        corrKurganovConvectionScheme(const corrKurganovConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const corrKurganovConvectionScheme&);

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > wIntermediate
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type, fvsPatchField, surfaceMesh>&
		) const;

		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > correction
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type, fvsPatchField, surfaceMesh>&
		) const;


public:

    //- Runtime type information
    TypeName("CorrKurganov");


    // Constructors

        //- Construct from flux and interpolation scheme
        corrKurganovConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type> >& scheme
        )
        :
            kurganovConvectionScheme<Type>(mesh, faceFlux, scheme)
        {}

        //- Construct from flux and Istream
        corrKurganovConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& flux,
            Istream& is
        )
        :
        	kurganovConvectionScheme<Type>(mesh,flux,is)
		{}


    // Member Functions

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

        /*tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > minMod
		(
			const GeometricField<Type, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type, fvsPatchField, surfaceMesh>&
		) const;*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "corrKurganovConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
