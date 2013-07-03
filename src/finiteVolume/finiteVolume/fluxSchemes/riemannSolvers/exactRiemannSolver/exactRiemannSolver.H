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
    Foam::fv::exactRiemannSolver

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    exactRiemannSolver.C

\*---------------------------------------------------------------------------*/

#ifndef exactRiemannSolver_H
#define exactRiemannSolver_H

#include <cmath>
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class exactRiemannSolver Declaration
\*---------------------------------------------------------------------------*/

class exactRiemannSolver

{
    // Private data

		const double& pL_;
		const double& pR_;
		const double& rhoL_;
		const double& rhoR_;
		const double& uL_;
		const double& uR_;
		const double gamma_;
		enum wave
		{
			NONE,
			SHOCK,
			RAREFACTION,
			VACUUM,
			CONTACT
		};
		exactRiemannSolver::wave waveL_;
		exactRiemannSolver::wave waveC_;
		exactRiemannSolver::wave waveR_;
		double pStar_;
		double uStarL_;
		double uStarR_;
		double rhoStarL_;
		double rhoStarR_;
		double speed_[5];

    // Private Member Functions

        //- Disallow default bitwise copy construct
        exactRiemannSolver(const exactRiemannSolver&);

        //- Disallow default bitwise assignment
        void operator=(const exactRiemannSolver&);

        double pStar
		(
			double
		) const;

		double f
		(
			const double&,
			const double&,
			const double&
		) const;

		double fDer
		(
			const double&,
			const double&,
			const double&
		) const;

		double uStar
		(
			double&
		) const;

		void uStar
		(
		);

		double pStarTwoRar
		(
		) const;

		void setWave
		(
			double&
		);

		void rhoStar
		(
		);

		void speed
		(
		);

		static Foam::word waveToString
		(
			const exactRiemannSolver::wave&
		);

		void check
		(
		)const;

protected:

public:

    //- Runtime type information
    TypeName("ExactReinmann");

    // Constructors

        //- Construct from flux and interpolation scheme
        exactRiemannSolver
        (
        	const double& pL,
        	const double& pR,
        	const double& rhoL,
        	const double& rhoR,
        	const double& uL,
        	const double& uR,
        	const double gamma
        )
    	:
    		pL_(pL),
    		pR_(pR),
    		rhoL_(rhoL),
    		rhoR_(rhoR),
    		uL_(uL),
    		uR_(uR),
    		gamma_(gamma),
			waveL_(NONE),
			waveC_(NONE),
			waveR_(NONE),
			pStar_(0.0),
			uStarL_(0.0),
			uStarR_(0.0),
			rhoStarL_(0.0),
			rhoStarR_(0.0)
        {}

        virtual ~exactRiemannSolver() {};
    // Public Member Functions

        static void test
        (
        );

        double p
		(
			double
		)const;

        double u
        (
        	double
        )const;

        double rho
        (
        	double
        )const;

        double speed
        (
        	int
        )const;

        void solve
		(
		);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



/*
#ifdef NoRepository
#   include "exactRiemannSolver.C"
#endif
*/



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
