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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
Foam::tmp<Foam::surfaceScalarField>
Foam::LimitedScheme<Type, Limiter, LimitFunc>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
	if(surfaceInterpolation::debug)
	{
		Info << "-----5----------" << endl;
	}
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                type() + "Limiter(" + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );
    surfaceScalarField& lim = tLimiter();

    tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh> >
        tlPhi = LimitFunc<Type>()(phi);

    const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
        lPhi = tlPhi();

    if(surfaceInterpolation::debug)
    {
    	Info << "-----5.1---------" << endl;
    }
    tmp<GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh> >
        tgradc(fvc::grad(lPhi));
    if(surfaceInterpolation::debug)
	{
		Info << "-----5.1---------" << endl;
	}

    const GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>&
        gradc = tgradc();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = lim.internalField();

    forAll(pLim, face)
    {
        label own = owner[face];
        label nei = neighbour[face];

        pLim[face] = Limiter::limiter
        (
            CDweights[face],
            this->faceFlux_[face],
            lPhi[own],
            lPhi[nei],
            gradc[own],
            gradc[nei],
            C[nei] - C[own]
        );
    }

    surfaceScalarField::GeometricBoundaryField& bLim = lim.boundaryField();

    forAll(bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            const Field<typename Limiter::phiType> plPhiP
            (
                lPhi.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::phiType> plPhiN
            (
                lPhi.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<typename Limiter::gradPhiType> pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );
            const Field<typename Limiter::gradPhiType> pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd = CDweights.boundaryField()[patchi].patch().delta();

            forAll(pLim, face)
            {
                pLim[face] = Limiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    plPhiP[face],
                    plPhiN[face],
                    pGradcP[face],
                    pGradcN[face],
                    pd[face]
                );
            }
        }
        else
        {
        	if (surfaceInterpolation::extrapolate)
        	{
        		if (surfaceInterpolation::debug)
        		{
        			Info << "extrapolating" << endl;
        		}
				// TODO verify if there is the possibility to use the coupled boundary instead of this implementation
				const scalarField& pCDweights = CDweights.boundaryField()[patchi];
				const scalarField& pFaceFlux =
					this->faceFlux_.boundaryField()[patchi];

				const Field<typename Limiter::phiType> plPhiP
				(
					lPhi.boundaryField()[patchi].patchInternalField()
				);
				const Field<typename Limiter::phiType> plPhiN
				(
					2*lPhi.boundaryField()[patchi] -
					lPhi.boundaryField()[patchi].patchInternalField()
				);
				const Field<typename Limiter::gradPhiType> pGradcP
				(
					gradc.boundaryField()[patchi].patchInternalField()
				);
				const Field<typename Limiter::gradPhiType> pGradcN
				(
					//TODO find a better model to guess the external gradient
					//-pGradcP		//specular trend on boundary
					pGradcP			//null second derivative on boundary
				);

				// Build the d-vectors
				vectorField pd = CDweights.boundaryField()[patchi].patch().delta();

				forAll(pLim, face)
				{
					pLim[face] = Limiter::limiter
					(
						pCDweights[face],
						pFaceFlux[face],
						plPhiP[face],
						plPhiN[face],
						pGradcP[face],
						pGradcN[face],
						pd[face]
					);
				}
        	}else
        	{
        		pLim = 1.0;
        	}
        }
    }

    return tLimiter;
}


// ************************************************************************* //
