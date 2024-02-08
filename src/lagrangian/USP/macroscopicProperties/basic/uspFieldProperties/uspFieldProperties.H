/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    uspFieldProperties

Description
    Stores all the information for the uspFieldProperties models
    defined within the fieldPropertiesDict, and selects & builds the models
    automatically.

\*---------------------------------------------------------------------------*/

#ifndef uspFieldProperties_H
#define uspFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "uspField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class uspFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class uspFieldProperties
{
    // Private data

        //- The entire field - properties dictionary
        //  (it contains multiple subDictionaries)

        const Time& time_;

        IOdictionary uspFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List<autoPtr<uspField>> fields_;


public:

    // Constructors

        //- Null Constructor
        uspFieldProperties
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for uspFoam
        uspFieldProperties
        (
            const Time& t,
            const polyMesh& mesh,
            uspCloud& cloud
        );

        //- Destructor
        ~uspFieldProperties() = default;


    // Member Functions

        void createFields();

        void calculateFields();

        void writeFields();

        inline List<autoPtr<uspField>>& fields();
        inline const List<autoPtr<uspField>>& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
