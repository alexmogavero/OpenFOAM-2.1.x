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
    Foam::fv::centralFluxScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    centralFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef centralConvectionScheme_H
#define centralConvectionScheme_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class centralFluxScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class centralFluxScheme
:
    public fv::fluxScheme<Type>
{
    // Private data

		const surfaceScalarField pos_;
		const surfaceScalarField neg_;
		Istream& is_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        centralFluxScheme(const centralFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const centralFluxScheme&);

protected:

        const basicThermo& thermo_;
		const volScalarField& rho_;
		const volVectorField& U_;


        virtual void weight
		(
		)=0;

		/*virtual void phiv
		(
		)=0;*/

		virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const=0;

		virtual tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const=0;

		virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > interpolate
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const=0;

		virtual tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > interpolate
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const=0;

		virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const=0;

		virtual tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const=0;

public:

    //- Runtime type information
    //TypeName("Central");

    // Constructors

        //- Construct from flux and interpolation scheme
    	centralFluxScheme
		(
			const fvMesh& mesh,
			const basicThermo& thermo,
			const volVectorField& U,
			const volScalarField& rho,
			Istream& is
		)
		:
			fluxScheme<Type>(mesh),
			pos_(IOobject("pos","dummy",mesh),mesh,dimensionedScalar("pos", dimless, 1.0)),
			neg_(IOobject("neg","dummy",mesh),mesh,dimensionedScalar("neg", dimless, -1.0)),
			is_(is),
			thermo_(thermo),
			rho_(rho),
			U_(U)
		{
		}

    // Member Functions

    	template<class Type1>
    	tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > interpolatePos
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

    	template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > interpolateNeg
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

    	void flux
		(
			surfaceScalarField& rhoFlux,
			surfaceVectorField& UFlux,
			surfaceScalarField& EFlux
		);

        /*tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
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
		) const;*/

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "centralFluxScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
