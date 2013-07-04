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
    Foam::fv::wasFluxScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    wasFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef wasFluxScheme_H
#define wasFluxScheme_H

//#include <cmath>
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "UList.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMesh.H"
#include "exactRiemannSolver.H"
#include "convectionScheme.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"
#include "fvcGrad.H"
#include "fluxScheme.H"
#include "wallFvPatch.H"
#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
#include "basicThermo.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class wasFluxScheme Declaration
\*---------------------------------------------------------------------------*/
template<class Type>
class wasFluxScheme
:
		public fv::fluxScheme<Type>
{
    // Private data
		const basicThermo& thermo_;
		const volScalarField& p_;
		volScalarField::GeometricBoundaryField pB_;
		volScalarField::GeometricBoundaryField pBB_;
		const volVectorField& U_;
		volVectorField::GeometricBoundaryField UB_;
		volVectorField::GeometricBoundaryField UBB_;
		const volScalarField& rho_;
		volScalarField::GeometricBoundaryField rhoB_;
		volScalarField::GeometricBoundaryField rhoBB_;
		surfaceScalarField p0_;
		surfaceScalarField::GeometricBoundaryField p0B_;
		surfaceScalarField p1_;
		surfaceScalarField::GeometricBoundaryField p1B_;
		surfaceScalarField p2_;
		surfaceScalarField::GeometricBoundaryField p2B_;
		surfaceScalarField p3_;
		surfaceScalarField::GeometricBoundaryField p3B_;
		surfaceScalarField rho0_;
		surfaceScalarField::GeometricBoundaryField rho0B_;
		surfaceScalarField rho1_;
		surfaceScalarField::GeometricBoundaryField rho1B_;
		surfaceScalarField rho2_;
		surfaceScalarField::GeometricBoundaryField rho2B_;
		surfaceScalarField rho3_;
		surfaceScalarField::GeometricBoundaryField rho3B_;
		surfaceScalarField u0_;
		surfaceScalarField::GeometricBoundaryField u0B_;
		surfaceScalarField u1_;
		surfaceScalarField::GeometricBoundaryField u1B_;
		surfaceScalarField u2_;
		surfaceScalarField::GeometricBoundaryField u2B_;
		surfaceScalarField u3_;
		surfaceScalarField::GeometricBoundaryField u3B_;
		surfaceScalarField speed0_;
		surfaceScalarField::GeometricBoundaryField speed0B_;
		surfaceScalarField speed1_;
		surfaceScalarField::GeometricBoundaryField speed1B_;
		surfaceScalarField speed2_;
		surfaceScalarField::GeometricBoundaryField speed2B_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        wasFluxScheme(const wasFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const wasFluxScheme&);

        void weight
        (
			double*,
			const double*,
			const double*,
			const double*,
			const double*,
			const vector&
        )const;

        void weight1
        (
        	double*,
        	const exactRiemannSolver&,
        	const exactRiemannSolver&,
        	const exactRiemannSolver&,
        	const vector&
        )const;

        double limiter
        (
        	const double r
        )const;

        double average
        (
        	const double*,
        	const double*
        )const;

        vector averageV
        (
        	const exactRiemannSolver&,
        	const vector&,
        	const vector&,
        	const tensor&,
        	const tensor&,
        	const vector&
        )const;

        vector averageV
        (
        	const scalar&,
        	const vector&,
        	const vector&,
        	const tensor&,
        	const tensor&,
        	const vector&
        )const;

        scalar rFunc
        (
			const scalar&,
			const scalar&,
			const scalar&,
			const scalar&
        ) const;

        void sampling
        (
        	double (exactRiemannSolver::*varFunc)(double)const,
        	const exactRiemannSolver&,
        	double*
        )const;

        void riemann
        (
        );

        void setBoundary
        (
			const volVectorField& gradP,
			const volTensorField& gradU,
			const volVectorField& gradRho
        );

protected:

public:

    //- Runtime type information
    TypeName("WAS");

    // Constructors

        wasFluxScheme
		(
			const fvMesh& mesh,
			const basicThermo& thermo,
			const volVectorField& U,
			const volScalarField& rho,
			Istream& schemeData
		)
		:
			fluxScheme<Type>(mesh),
			thermo_(thermo),
			p_(thermo.p()),
			pB_(p_.boundaryField()),
			pBB_(p_.boundaryField()),
			U_(U),
			UB_(U_.boundaryField()),
			UBB_(U_.boundaryField()),
			rho_(rho),
			rhoB_(rho_.boundaryField()),
			rhoBB_(rho_.boundaryField()),
			p0_
			(
				IOobject
				(
					"p1",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar("p1", p_.dimensions(), 0.0)
			),
			p0B_(p0_.boundaryField()),
			p1_(p0_),
			p1B_(p1_.boundaryField()),
			p2_(p0_),
			p2B_(p2_.boundaryField()),
			p3_(p0_),
			p3B_(p3_.boundaryField()),
			rho0_
			(
				IOobject
				(
					"rho1",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar("rho1", rho_.dimensions(), 0.0)
			),
			rho0B_(rho0_.boundaryField()),
			rho1_(rho0_),
			rho1B_(rho1_.boundaryField()),
			rho2_(rho0_),
			rho2B_(rho2_.boundaryField()),
			rho3_(rho0_),
			rho3B_(rho3_.boundaryField()),
			u0_
			(
				IOobject
				(
					"u1",
					mesh.time().timeName(),
					mesh
				),
				mesh,
				dimensionedScalar("u1", U_.dimensions(), 0.0)
			),
			u0B_(u0_.boundaryField()),
			u1_(u0_),
			u1B_(u1_.boundaryField()),
			u2_(u0_),
			u2B_(u2_.boundaryField()),
			u3_(u0_),
			u3B_(u3_.boundaryField()),
			speed0_(u0_),
			speed0B_(speed0_.boundaryField()),
			speed1_(u0_),
			speed1B_(speed1_.boundaryField()),
			speed2_(u0_),
			speed2B_(speed2_.boundaryField())
		{}

        virtual ~wasFluxScheme();
    // Public Member Functions
        void test
		(
		);

        void calculate
		(
			surfaceScalarField&,
			surfaceScalarField&,
			surfaceVectorField&
		);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#ifdef NoRepository
#   include "wasFluxScheme.C"
#endif



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
