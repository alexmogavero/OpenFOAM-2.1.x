/*
 * cellDebug.H
 *
 *  Created on: 16 Apr 2013
 *      Author: trb12187
 */


const labelList& faceOwner = mesh.faceOwner();
const labelList& faceNeigh = mesh.neighbour();
const polyBoundaryMesh& boundary = mesh.boundaryMesh();
const int cellId = 119;

const surfaceScalarField& F1 = phiEp;
const surfaceScalarField& F2 = phi;

Info << "Cell faces detail for cell# " << cellId << endl;
Info << "ID\tinfo\t\t\tphiEp1\t\t\tphiEp2\t\t\trho_pos\t\trho_neg\t\te_pos\t\te_neg" << endl;

forAll(faceOwner,i)
{
	if(faceOwner[i]==cellId)
	{
		if(mesh.isInternalFace(i))
		{
			Info << i << "\t" << "Owner " << faceNeigh[i] << "\t";
			Info << "\t" << F1[i];
			Info << "\t" << F2[i];
			Info << "\t" << rho_pos[i];
			Info << "\t" << rho_neg[i];
			Info << "\t" << e_pos[i];
			Info << "\t" << e_neg[i];
			Info << endl;
		}else
		{
			forAll(boundary[boundary.whichPatch(i)],ii)
			{
				if(boundary[boundary.whichPatch(i)].faceCells()[ii]==cellId)
				{
					Info << ii << "\t" << boundary.names()[boundary.whichPatch(i)];
					Info << "\t" << F1.boundaryField()[boundary.whichPatch(i)][ii];
					Info << "\t" << F2.boundaryField()[boundary.whichPatch(i)][ii];
					Info << "\t" << rho_pos.boundaryField()[boundary.whichPatch(i)][ii];
					Info << "\t" << rho_neg.boundaryField()[boundary.whichPatch(i)][ii];
					Info << "\t" << e_pos.boundaryField()[boundary.whichPatch(i)][ii];
					Info << "\t" << e_neg.boundaryField()[boundary.whichPatch(i)][ii];
					Info << endl;
				}
			}
		}
	}
}

forAll(faceNeigh,i)
{
	if(faceNeigh[i]==cellId)
	{
		if(mesh.isInternalFace(i))
		{
			Info << i << "\t" << "Neighbour " << faceOwner[i] << "\t";
			Info << "\t" << -F1[i];
			Info << "\t" << -F2[i];
			Info << "\t" << rho_pos[i];
			Info << "\t" << rho_neg[i];
			Info << "\t" << e_pos[i];
			Info << "\t" << e_neg[i];
			Info << endl;
		}
	}
}
//volScalarField divphiEp_e(fvc::div(phiEp_e));
DimensionedField<double,volMesh> volume(mesh.V());
Info << "divphiEp " << divphiEp[cellId]*volume[cellId] << " divphiEp2 " << divphiEp2[cellId]*volume[cellId] << endl;

//face interpolation comparison
const surfaceScalarField& var_pos = e_pos;
const surfaceScalarField& var_neg = e_neg;
const volScalarField& var = e;
const volVectorField varG = fvc::grad(e);
//surfaceInterpolation::extrapolate = false;
tmp<limitedSurfaceInterpolationScheme<scalar> > scheme = limitedSurfaceInterpolationScheme<scalar>::
		New(mesh,pos,mesh.interpolationScheme("reconstruct(T)"));
scheme().extrapolate_ = true;
const tmp<surfaceScalarField> l = scheme().limiter(var);
//surfaceInterpolation::extrapolate = false;
const int faceID [3] = {237,235,233};
const int cellID [4] = {119,118,117,116};

Info << "B   " << mesh.Cf().boundaryField()[mesh.boundary().findPatchID("wall")][0][0] << " "
     << var.boundaryField()[mesh.boundary().findPatchID("wall")][0] << endl;

Info << "0   " << mesh.Cf().boundaryField()[mesh.boundary().findPatchID("wall")][0][0] << " "
     << var_pos.boundaryField()[mesh.boundary().findPatchID("wall")][0] << " ";
Info << var_neg.boundaryField()[mesh.boundary().findPatchID("wall")][0] << " ";
Info << var_pos.boundaryField()[mesh.boundary().findPatchID("wall")][0] -
		var_neg.boundaryField()[mesh.boundary().findPatchID("wall")][0] << " "
     << l().boundaryField()[mesh.boundary().findPatchID("wall")][0] << endl;

Info << cellID[0] << " " << mesh.C()[cellID[0]][0] << " " << var[cellID[0]]
     << " " << varG[cellID[0]][0] << endl;

int i;
for(i=0;i<3;i++)
{
	Info << faceID[i] << " " << mesh.Cf()[faceID[i]][0] << " " << var_pos[faceID[i]] << " " << var_neg[faceID[i]] <<
			" " << var_pos[faceID[i]] - var_neg[faceID[i]] <<
			" " << l()[faceID[i]] << endl;
	Info << cellID[i+1] << " " << mesh.C()[cellID[i+1]][0] << " " << var[cellID[i+1]] << " " << varG[cellID[i+1]][0] << endl;
}
