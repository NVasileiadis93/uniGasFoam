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
    Foam::variableHardSphere

Description
    This class provides functions to use the variable hard sphere binary collision model.

    Use in the uniGasProperties file.

    Reference:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (pages 41, 97, 222, 426, 428).
        
Usage
    Example of the binary collision model specification:

    \verbatim
    binaryCollisionModel            variableHardSphere;
    collisionProperties
    {
        Tref                        273;
    }
    \endverbatim

SourceFiles
    variableHardSphere.C

\*---------------------------------------------------------------------------*/

#ifndef variableHardSphere_H
#define variableHardSphere_H

#include "binaryCollisionModel.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class variableHardSphere Declaration
\*---------------------------------------------------------------------------*/

class variableHardSphere
:
    public binaryCollisionModel
{
    // Private data

        //- The properties dictionary
        const dictionary propertiesDict_;

        //- Reference temperature
        const scalar Tref_;


public:

    //- Runtime type information
    TypeName("variableHardSphere");

    //- Construct from dictionary
    variableHardSphere
    (
        const dictionary& dict,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~variableHardSphere() = default;


    // Member Functions

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

        virtual const dictionary& propertiesDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "variableHardSphere.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
