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
    Foam::fv::kurganovFluxScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    kurganovFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef kurganovFluxScheme_H
#define kurganovFluxScheme_H

#include "centralFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class kurganovFluxScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class kurganovFluxScheme
:
    public fv::centralFluxScheme<Type>
{
    // Private data

		/*const basicThermo& thermo_;
		const volScalarField& rho_;
		const volVectorField& U_;*/


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kurganovFluxScheme(const kurganovFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const kurganovFluxScheme&);

protected:

        surfaceScalarField ap_;
		surfaceScalarField am_;
		surfaceScalarField phiv_pos_;
		surfaceScalarField phiv_neg_;
		surfaceScalarField a_pos_;
		surfaceScalarField a_neg_;
		surfaceScalarField aphiv_pos_;
		surfaceScalarField aphiv_neg_;
		surfaceScalarField aSf_;

    // Protected member functions
        void weight
		(
		);

        template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > fluxT
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

		virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const;

		virtual tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > flux
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const;

        template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > interpolateT
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > interpolate
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > interpolate
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const;

        template<class Type1>
		tmp<GeometricField<Type1, fvsPatchField, surfaceMesh> > diffusionT
		(
			const GeometricField<Type1, fvPatchField, volMesh>&
		) const;

		virtual tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<scalar, fvPatchField, volMesh>&
		) const;

		virtual tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > diffusion
		(
			const GeometricField<vector, fvPatchField, volMesh>&
		) const;

public:

    //- Runtime type information
    TypeName("Kurganov");


    // Constructors

        //- Construct
        kurganovFluxScheme
		(
			const fvMesh& mesh,
			const basicThermo& thermo,
			const volVectorField& U,
			const volScalarField& rho,
			Istream& is
		)
        :
        	centralFluxScheme<Type>(mesh,thermo,U,rho,is),
        	/*thermo_(thermo),
        	rho_(rho),
        	U_(U),*/
        	ap_(IOobject("ap","dummy",mesh),mesh,dimensionedScalar("ap", dimVolume/dimTime, 0.0)),
        	am_(ap_),
        	phiv_pos_(ap_),
        	phiv_neg_(ap_),
        	a_pos_(IOobject("ap","dummy",mesh),mesh,dimensionedScalar("ap", dimless, 0.0)),
        	a_neg_(a_pos_),
        	aphiv_pos_(ap_),
        	aphiv_neg_(ap_),
        	aSf_(ap_)
		{}


    // Member Functions

        /*void flux
        (
        	surfaceScalarField& rhoFlux,
        	surfaceVectorField& UFlux,
        	surfaceScalarField& EFlux
        );*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kurganovFluxScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
