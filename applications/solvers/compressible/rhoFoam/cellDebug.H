/*
 * cellDebug.H
 *
 *  Created on: 16 Apr 2013
 *      Author: trb12187
 */


Info << "Proc: " << Pstream::myProcNo() << endl;
const labelList& faceOwner = mesh.faceOwner();
const labelList& faceNeigh = mesh.neighbour();
const polyBoundaryMesh& boundary = mesh.boundaryMesh();
const int cellId = 92377; //119; //38336;

const surfaceScalarField& F1 = rhoFlux;
const surfaceScalarField F2 = UFlux.component(vector::X);
const surfaceScalarField& Q1 = rhoAve;
const surfaceScalarField Q2 = UAve & (mesh.Sf()/mesh.magSf());
const surfaceScalarField Q3 = mesh.magSf();
Info << "Cell faces detail for cell# " << cellId << endl;
Info << "ID\tinfo\t\t\tF1\t\t\tF2\t\t\tQ1\t\tQ2\t\tQ3 " << endl;

forAll(faceOwner,i)
{
	if(faceOwner[i]==cellId)
	{
		if(mesh.isInternalFace(i))
		{
			Info << i << "\t" << "Owner " << faceNeigh[i] << "\t";
			Info << "\t" << F1[i];
			Info << "\t" << F2[i];
			Info << "\t" << Q1[i];
			Info << "\t" << Q2[i];
			Info << "\t" << Q3[i];
			//Info << "\t" << e_neg[i];
			Info << endl;
		}else
		{
			forAll(boundary[boundary.whichPatch(i)],ii)
			{
				if(boundary[boundary.whichPatch(i)].faceCells()[ii]==cellId)
				{
					if(F1.boundaryField()[boundary.whichPatch(i)].size()>0)
					{
						Info << ii << "\t" << boundary.names()[boundary.whichPatch(i)];
						Info << "\t" << F1.boundaryField()[boundary.whichPatch(i)][ii];
						Info << "\t" << F2.boundaryField()[boundary.whichPatch(i)][ii];
						Info << "\t" << Q1.boundaryField()[boundary.whichPatch(i)][ii];
						Info << "\t" << Q2.boundaryField()[boundary.whichPatch(i)][ii];
						Info << "\t" << Q3.boundaryField()[boundary.whichPatch(i)][ii];
						//Info << "\t" << e_neg.boundaryField()[boundary.whichPatch(i)][ii];
						Info << endl;
					}
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
			Info << "\t" << Q1[i];
			Info << "\t" << Q2[i];
			Info << "\t" << Q3[i];
			//Info << "\t" << e_neg[i];
			Info << endl;
		}
	}
}
volScalarField divF1(fvc::div(F1));
volScalarField divF2(fvc::div(F2));
DimensionedField<double,volMesh> volume(mesh.V());
Info << "divF1 " << divF1[cellId]*volume[cellId] << " divF2 " << divF2[cellId]*volume[cellId] << endl;

