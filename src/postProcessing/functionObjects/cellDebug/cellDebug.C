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

#include "cellDebug.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::cellDebug, 0);

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::cellDebug::modeType,
        2
    >::names[] =
    {
        "magnitude",
        "component"
    };
}


const Foam::NamedEnum<Foam::cellDebug::modeType, 2>
Foam::cellDebug::modeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellDebug::cellDebug
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    write_(true),
    log_(false),
    mode_(mdMag),
    fieldSet_(),
    cellDebugFilePtr_(NULL),
    cellId_(0),
    procNo_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "cellDebug::cellDebug"
            "(const objectRegistry& obr, const dictionary& dict)"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellDebug::~cellDebug()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellDebug::read(const dictionary& dict)
{
    if (active_)
    {
        write_ = dict.lookupOrDefault<Switch>("write", true);
        log_ = dict.lookupOrDefault<Switch>("log", false);

        mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "magnitude")];
        dict.lookup("fields") >> fieldSet_;
        cellId_ = dict.lookupOrDefault<label>("cell",0);
        procNo_ = dict.lookupOrDefault<label>("processor",0);
    }
}


void Foam::cellDebug::makeFile()
{
    // Create the cellDebug file if not already created
    if (cellDebugFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating cellDebug file." << endl;
        }

        // File update
        if (Pstream::myProcNo()==procNo_) //(Pstream::master())
        {
            fileName cellDebugDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                cellDebugDir =
                    obr_.time().path()/".."/name_/obr_.time().timeName();
            }
            else
            {
                cellDebugDir =
                    obr_.time().path()/name_/obr_.time().timeName();
            }

            // Create directory if does not exist.
            mkDir(cellDebugDir);

            // Open new file at start up
            cellDebugFilePtr_.reset
            (
                new OFstream(cellDebugDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::cellDebug::writeFileHeader()
{
    if (cellDebugFilePtr_.valid())
    {
    	cellDebugFilePtr_() << "Cell faces detail for cell# " << cellId_
    					<< " within proc domain# " << Pstream::myProcNo() << endl;
        /*cellDebugFilePtr_()
            << "# Time" << tab << "field" << tab << "min" << tab << "max"
            << endl;*/
    }
}


void Foam::cellDebug::execute()
{
    // Do nothing - only valid on write
}


void Foam::cellDebug::end()
{
    // Do nothing - only valid on write
}


void Foam::cellDebug::write()
{
    if (active_ & write_)
    {
        // Create the cellDebug file if not already created
        makeFile();

		const surfaceScalarField& F1 = obr_.lookupObject<surfaceScalarField>(fieldSet_[0]);
		const fvMesh& mesh = F1.mesh();
		const labelList& faceOwner = mesh.faceOwner();
		const labelList& faceNeigh = mesh.neighbour();
		const polyBoundaryMesh& boundary = mesh.boundaryMesh();

		if (Pstream::myProcNo()==procNo_)
		{
			cellDebugFilePtr_() << "Time= " << obr_.time().value() << endl;
			cellDebugFilePtr_() << "face" << tab << "cells";
			forAll(fieldSet_, fieldI)
			{
				const surfaceScalarField& field = obr_.lookupObject<surfaceScalarField>(fieldSet_[fieldI]);
				cellDebugFilePtr_() << tab << field.name();
			}
			cellDebugFilePtr_() << endl;

			forAll(faceOwner,i)
			{
				if(faceOwner[i]==cellId_)
				{
					if(mesh.isInternalFace(i))
					{
						cellDebugFilePtr_() << i << tab << cellId_ << "->" << faceNeigh[i];
						forAll(fieldSet_, fieldI)
						{
							const surfaceScalarField& field = obr_.lookupObject<surfaceScalarField>(fieldSet_[fieldI]);
							cellDebugFilePtr_() << tab << field[i];
						}
						cellDebugFilePtr_() << endl;
					}else
					{
						forAll(boundary[boundary.whichPatch(i)],ii)
						{
							if(boundary[boundary.whichPatch(i)].faceCells()[ii]==cellId_)
							{
								if(F1.boundaryField()[boundary.whichPatch(i)].size()>0)
								{
									cellDebugFilePtr_() << ii << tab << boundary.names()[boundary.whichPatch(i)];
									forAll(fieldSet_, fieldI)
									{
										const surfaceScalarField& field = obr_.lookupObject<surfaceScalarField>(fieldSet_[fieldI]);
										cellDebugFilePtr_() << tab << field.boundaryField()[boundary.whichPatch(i)][ii];
									}
									cellDebugFilePtr_() << endl;
								}
							}
						}
					}
				}
			}

			forAll(faceNeigh,i)
			{
				if(faceNeigh[i]==cellId_)
				{
					if(mesh.isInternalFace(i))
					{
						cellDebugFilePtr_() << i << tab << faceOwner[i] << "->" << cellId_;
						forAll(fieldSet_, fieldI)
						{
							const surfaceScalarField& field = obr_.lookupObject<surfaceScalarField>(fieldSet_[fieldI]);
							cellDebugFilePtr_() << tab << field[i];
						}
						cellDebugFilePtr_() << endl;
					}
				}
			}
			cellDebugFilePtr_() << endl;
		}
    }
}


template<>
void Foam::cellDebug::calcMinMaxFields<Foam::scalar>
(
    const word& fieldName
)
{
    if (obr_.foundObject<volScalarField>(fieldName))
    {
        const volScalarField& field =
            obr_.lookupObject<volScalarField>(fieldName);
        const scalar minValue = min(field).value();
        const scalar maxValue = max(field).value();

        if (Pstream::master())
        {
            if (write_)
            {
                cellDebugFilePtr_()
                    << obr_.time().value() << tab
                    << fieldName << tab << minValue << tab << maxValue << endl;
            }

            if (log_)
            {
                Info<< "cellDebug output:" << nl
                    << "    min(" << fieldName << ") = " << minValue << nl
                    << "    max(" << fieldName << ") = " << maxValue << nl
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
