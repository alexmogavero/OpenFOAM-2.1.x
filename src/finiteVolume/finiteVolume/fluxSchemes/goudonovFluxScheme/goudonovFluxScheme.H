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
    Foam::fv::goudonovFluxScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    goudonovFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef goudonovFluxScheme_H
#define goudonovFluxScheme_H

//#include <cmath>
#include "typeInfo.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMesh.H"
#include "exactRiemannSolver.H"
#include "fvcGrad.H"
#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class goudonovFluxScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class goudonovFluxScheme
:
		public fv::fluxScheme<Type> //refCount
{
    // Private data
		const volScalarField& p_;
		const volVectorField& U_;
		const volScalarField& rho_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        goudonovFluxScheme(const goudonovFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const goudonovFluxScheme&);

protected:

public:

    //- Runtime type information
    TypeName("Goudonov");

    // Constructors

        goudonovFluxScheme
        (
        	const fvMesh& mesh,
        	//const surfaceScalarField& flux,
        	const volScalarField& p,
        	const volVectorField& U,
        	const volScalarField& rho,
        	Istream& schemeData
        )
		:
			fluxScheme<Type>(mesh),
			p_(p),
			U_(U),
			rho_(rho)
        {}

        virtual ~goudonovFluxScheme();

    // Public Member Functions

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
#   include "goudonovFluxScheme.C"
#endif



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
