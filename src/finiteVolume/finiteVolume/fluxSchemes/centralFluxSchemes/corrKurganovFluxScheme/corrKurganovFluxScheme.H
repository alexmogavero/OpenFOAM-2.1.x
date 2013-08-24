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
    Foam::fv::corrKurganovFluxScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    corrKurganovFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef corrKurganovFluxScheme_H
#define corrKurganovFluxScheme_H

#include "kurganovFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class corrKurganovFluxScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class corrKurganovFluxScheme
:
    public fv::kurganovFluxScheme<Type>
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        corrKurganovFluxScheme(const corrKurganovFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const corrKurganovFluxScheme&);

        template<class Type1>
        tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > wIntermediate
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type1, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type1, fvsPatchField, surfaceMesh>&
		) const;

        template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > correction
		(
			const surfaceScalarField&,
			const surfaceScalarField&,
			const GeometricField<Type1, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type1, fvsPatchField, surfaceMesh>&
		) const;

		template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > fluxT
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const;

		template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > diffusionT
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const;


public:

    //- Runtime type information
    TypeName("CorrKurganov");


    // Constructors

        //- Construct from flux and interpolation scheme
    	corrKurganovFluxScheme
		(
			const fvMesh& mesh,
			const basicThermo& thermo,
			const volVectorField& U,
			const volScalarField& rho,
			Istream& is
		)
		:
			kurganovFluxScheme<Type>(mesh,thermo,U,rho,is)
        {}


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "corrKurganovFluxScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
