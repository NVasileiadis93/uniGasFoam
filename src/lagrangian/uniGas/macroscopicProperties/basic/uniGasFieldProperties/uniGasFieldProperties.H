/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::uniGasFieldProperties

Description
    Stores all the information for the uniGasFieldProperties models
    defined within the fieldPropertiesDict, and selects & builds the models
    automatically.

Source Files
    uniGasFieldPropertiesI.H
    uniGasFieldProperties.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasFieldProperties_H
#define uniGasFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "uniGasField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class uniGasFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class uniGasFieldProperties
{
    // Private data

        //- The entire field - properties dictionary
        //  (it contains multiple subDictionaries)

        const Time& time_;

        IOdictionary uniGasFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List<autoPtr<uniGasField>> fields_;


public:

    // Constructors

        //- Null Constructor
        uniGasFieldProperties
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for uniGasFoam
        uniGasFieldProperties
        (
            const Time& t,
            const polyMesh& mesh,
            uniGasCloud& cloud
        );

        //- Destructor
        ~uniGasFieldProperties() = default;


    // Member Functions

        void createFields();

        void calculateFields();

        void writeFields();

        inline List<autoPtr<uniGasField>>& fields();
        inline const List<autoPtr<uniGasField>>& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniGasFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
