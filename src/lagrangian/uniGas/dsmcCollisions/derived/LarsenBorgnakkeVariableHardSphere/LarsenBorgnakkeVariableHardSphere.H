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
    Foam::LarsenBorgnakkeVariableHardSphere

Description
    Variable Hard Sphere dsmc collision model with Larsen Borgnakke internal
    energy redistribution.

    Use in the uniGasProperties file.

Usage
    Example of the dsmc collision model specification:

    \verbatim
    dsmcCollisionModel        LarsenBorgnakkeVariableHardSphere;
    collisionProperties
    {
        Tref                                273;
        rotationalRelaxationCollisionNumber 5;
        electronicRelaxationCollisionNumber 50;
    }
    \endverbatim

SourceFiles
    LarsenBorgnakkeVariableHardSphere.C

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVariableHardSphere_H
#define LarsenBorgnakkeVariableHardSphere_H

#include "dsmcCollisionModel.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class LarsenBorgnakkeVariableHardSphere Declaration
\*---------------------------------------------------------------------------*/

class LarsenBorgnakkeVariableHardSphere
:
    public dsmcCollisionModel
{
    // Private data

        //- The properties dictionary
        const dictionary propertiesDict_;

        //- Reference temperature
        const scalar Tref_;

        //- Rotational relaxation collision number
        const scalar rotationalRelaxationCollisionNumber_;

        //- Electronic relaxation collision number
        const scalar electronicRelaxationCollisionNumber_;


public:

    //- Runtime type information
    TypeName("LarsenBorgnakkeVariableHardSphere");

    //- Construct from dictionary
    LarsenBorgnakkeVariableHardSphere
    (
        const dictionary& dict,
        uniGasCloud& cloud
    );


    // Destructor
    virtual ~LarsenBorgnakkeVariableHardSphere() = default;


    // Member Functions

        //- Return the coefficients dictionary
        virtual const dictionary& propertiesDict() const;

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Return the collision cross section * relative velocity product
        virtual scalar sigmaTcR
        (
            const uniGasParcel& pP,
            const uniGasParcel& pQ
        ) const;

        //- Apply collision
        virtual void collide
        (
            uniGasParcel& pP,
            uniGasParcel& pQ,
            label& cellI
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
