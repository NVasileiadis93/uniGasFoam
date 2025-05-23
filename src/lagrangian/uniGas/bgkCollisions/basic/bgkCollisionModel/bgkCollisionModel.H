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
    Foam::bgkCollisionModel

Description
    Templated particle bgk collision class

SourceFiles
    bgkCollisionModel.C

\*---------------------------------------------------------------------------*/

#ifndef bgkCollisionModel_H
#define bgkCollisionModel_H

#include "fvMesh.H"
#include "uniGasParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uniGasCloud;
class uniGasParcel;

/*---------------------------------------------------------------------------*\
                    Class bgkCollisionModel Declaration
\*---------------------------------------------------------------------------*/

class bgkCollisionModel
{
    // Private data

protected:

        //- The cloud dictionary
        const dictionary& dict_;

        //- reference to the mesh
        const fvMesh& mesh_;

        // reference to the owner cloud class
        uniGasCloud& cloud_;


public:

    //- Runtime type information
    TypeName("bgkCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        bgkCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uniGasCloud& owner
        ),
        (dict, mesh, owner)
    );


    // Constructors

        //- Construct from components
        bgkCollisionModel
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uniGasCloud& owner
        );


    //- Destructor
    virtual ~bgkCollisionModel() = default;


    //- Selector
    static autoPtr<bgkCollisionModel> New
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& owner
    );


    // Access

        //- Return the dictionary
        const dictionary& dict() const;


    // Member Functions

        virtual void collide() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
