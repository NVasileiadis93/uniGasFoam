/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::relaxationCollisionModel

Description
    Templated particle relaxation class

SourceFiles
    relaxationCollisionModel.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationCollisionModel_H
#define relaxationCollisionModel_H

#include "fvMesh.H"
#include "uspParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uspCloud;
class uspParcel;

/*---------------------------------------------------------------------------*\
                    Class relaxationCollisionModel Declaration
\*---------------------------------------------------------------------------*/

class relaxationCollisionModel
{
    // Private data

protected:

        //- The cloud dictionary
        const dictionary& dict_;

        //- reference to the mesh
        const fvMesh& mesh_;

        // reference to the owner cloud class
        uspCloud& cloud_;


public:

    //- Runtime type information
    TypeName("relaxationCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relaxationCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uspCloud& owner
        ),
        (dict, mesh, owner)
    );


    // Constructors

        //- Construct from components
        relaxationCollisionModel
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uspCloud& owner
        );


    //- Destructor
    virtual ~relaxationCollisionModel() = default;


    //- Selector
    static autoPtr<relaxationCollisionModel> New
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& owner
    );


    // Access

        //- Return the dictionary
        const dictionary& dict() const;


    // Member Functions

        virtual void relax() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
