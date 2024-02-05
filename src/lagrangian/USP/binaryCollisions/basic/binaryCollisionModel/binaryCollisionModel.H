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
    Foam::binaryCollisionModel

Description
    Templated USP particle collision class

SourceFiles
    binaryCollisionModel.C
    binaryCollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef binaryCollisionModel_H
#define binaryCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uspCloud;
class uspParcel;

/*---------------------------------------------------------------------------*\
                    Class binaryCollisionModel Declaration
\*---------------------------------------------------------------------------*/

class binaryCollisionModel
{
    // Private data

protected:

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        uspCloud& cloud_;


public:

    //- Runtime type information
    TypeName("binaryCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        binaryCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            uspCloud& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        binaryCollisionModel(uspCloud& owner);

        //- Construct from components
        binaryCollisionModel
        (
            const dictionary& dict,
            uspCloud& owner
        );


    //- Destructor
    virtual ~binaryCollisionModel() = default;


    //- Selector
    static autoPtr<binaryCollisionModel> New
    (
        const dictionary& dict,
        uspCloud& owner
    );


    // Access

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        virtual const dictionary& propertiesDict() const = 0;


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const = 0;

        //- Return the collision cross section * relative velocity product
        virtual scalar sigmaTcR
        (
            const uspParcel& pP,
            const uspParcel& pQ
        ) const = 0;

        //- Apply collision
        virtual void collide
        (
            uspParcel& pP,
            uspParcel& pQ,
            label& cellI
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
