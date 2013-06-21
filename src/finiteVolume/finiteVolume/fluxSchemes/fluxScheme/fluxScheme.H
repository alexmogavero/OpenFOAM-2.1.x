/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute itand/or modify it
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
    Foam::fv::fluxScheme

Description
    Abstract base class for convection schemes.

SourceFiles
    fluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_H
#define fluxScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "multivariateSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxScheme
:
    public refCount
{
    // Private data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fluxScheme,
            Istream,
            (
                const fvMesh& mesh,
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& rho,
                Istream& schemeData
            ),
            (mesh, p, U, rho, schemeData)
        );

        /*declareRunTimeSelectionTable
        (
            tmp,
            fluxScheme,
            Multivariate,
            (
                const fvMesh& mesh,
                const typename multivariateSurfaceInterpolationScheme<Type>::
                    fieldTable& fields,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, fields, faceFlux, schemeData)
        );*/


    // Constructors

        //- Copy construct
        fluxScheme(const fluxScheme&);

        //- Construct from mesh, flux and Istream
        fluxScheme
        (
            const fvMesh& mesh
        )
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new fluxScheme created on freestore
        static tmp<fluxScheme<Type> > New
        (
            const fvMesh& mesh,
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& rho,
            Istream& schemeData
        );


        //- Return a pointer to a new multivariate fluxScheme
        //  created on freestore
        /*static tmp<fluxScheme > New
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );*/


    //- Destructor
    virtual ~fluxScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual void
        calculate
		(
			surfaceScalarField&,
			surfaceScalarField&,
			surfaceVectorField&
		) = 0;

        /*virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;*/

        /*virtual tmp<fvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;*/


    // Member operators

        void operator=(const fluxScheme<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvFluxTypeScheme(SS,Type)                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
fluxScheme<Type>::addIstreamConstructorToTable<SS<Type> >                      \
    add##SS##Type##IstreamConstructorToTable_;

/*
#define makeFvFluxScheme(SS)                                                   \
                                                                               \
defineTypeNameAndDebug(SS, 0);                                                 \
                                                                               \
fluxScheme::addIstreamConstructorToTable<SS >                                  \
    add##SS##IstreamConstructorToTable_;
*/

#define makeFvFluxScheme(SS)                                                   \
makeFvFluxTypeScheme(SS, scalar)                                               \
makeFvFluxTypeScheme(SS, vector)                                       \
makeFvFluxTypeScheme(SS, sphericalTensor)                                \
makeFvFluxTypeScheme(SS, symmTensor)                                     \
makeFvFluxTypeScheme(SS, tensor)


/*
#define makeMultivariateFvConvectionTypeScheme(SS, Type)                       \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
fluxScheme<Type>::addMultivariateConstructorToTable<SS<Type> >           \
    add##SS##Type##MultivariateConstructorToTable_;


#define makeMultivariateFvfluxScheme(SS)                                 \
                                                                               \
makeMultivariateFvConvectionTypeScheme(SS, scalar)                             \
makeMultivariateFvConvectionTypeScheme(SS, vector)                             \
makeMultivariateFvConvectionTypeScheme(SS, sphericalTensor)                    \
makeMultivariateFvConvectionTypeScheme(SS, symmTensor)                         \
makeMultivariateFvConvectionTypeScheme(SS, tensor)*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


   //Uncomment if the class is template


#ifdef NoRepository
#   include "fluxScheme.C"
#endif





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
