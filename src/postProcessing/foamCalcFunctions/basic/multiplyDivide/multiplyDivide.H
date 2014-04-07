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
    Foam::calcTypes::multiplyDivide

Description
    adds/subtracts a field or value to/from a base field.

    New field name specified by -resultName option, or automatically as:
        \<baseFieldName\>_add_<multiplyDivideFieldName>
        \<baseFieldName\>_add_value
        \<baseFieldName\>_subtract_<multiplyDivideFieldName>
        \<baseFieldName\>_subtract_value

    Example usage:
        multiplyDivide p add -value 100000 -resultName pAbs
        multiplyDivide U subtract -field U0

SourceFiles
    multiplyDivide.C
    writemultiplyDivideField.C
    writemultiplyDivideValue.C

\*---------------------------------------------------------------------------*/

#ifndef multiplyDivide_H
#define multiplyDivide_H

#include "calcType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                          Class multiplyDivide Declaration
\*---------------------------------------------------------------------------*/

class multiplyDivide
:
    public calcType
{
public:

    enum calcTypes
    {
        FIELD,
        VALUE
    };

    enum calcModes
    {
        MULTIPLY,
        DIVIDE
    };


private:

    // Private data

        //- Name of base field (to multiplyDivide to)
        word baseFieldName_;

        //- Calc type as given by enumerations above
        calcTypes calcType_;

        //- Name of field to add/subtract
        word multiplyDivideFieldName_;

        //- String representation of value to add/subtract
        string multiplyDivideValueStr_;

        //- Name of result field
        word resultName_;

        //- Mode - multiplyDivide/subtract
        calcModes calcMode_;


    // Private Member Functions

        // Output

            //- Calc and output field multiplyDivideitions
            void writeMultiplyDivideFields
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );

            //- Calc and output field and value multiplyDivideitions
            void writeMultiplyDivideValues
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );


        //- Disallow default bitwise copy construct
        multiplyDivide(const multiplyDivide&);

        //- Disallow default bitwise assignment
        void operator=(const multiplyDivide&);


protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );


        // I-O

            //- Write multiplyDivide field
            template<class Type>
            void writeMultiplyDivideField
            (
                const IOobject& baseHeader,
                const IOobject& multiplyDivideHeader,
                const fvMesh& mesh,
                bool& processed
            );

            //- Write multiplyDivide value
            template<class Type>
            void writeMultiplyDivideValue
            (
                const IOobject& baseHeader,
                const string& valueStr,
                const fvMesh& mesh,
                bool& processed
            );


public:

    //- Runtime type information
    TypeName("multiplyDivide");


    // Constructors

        //- Construct null
        multiplyDivide();


    //- Destructor
    virtual ~multiplyDivide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "writeMultiplyDivideField.C"
#   include "writeMultiplyDivideValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
