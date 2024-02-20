/*---------------------------------------------------------------------------* \
  =========                 |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::variableSoftSphere

Description
    This class provides functions to use the variable soft sphere binary collision model.

    Use in the uspProperties file.

    Reference:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (pages 41, 97, 222, 426, 428).

Usage
    Example of the binary collision model specification:

    \verbatim
    binaryCollisionModel            variableSoftSphere;
    collisionProperties
    {
        Tref                        273;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef variableSoftSphere_H
#define variableSoftSphere_H

#include "binaryCollisionModel.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class variableSoftSphere Declaration
\*---------------------------------------------------------------------------*/

class variableSoftSphere
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
    TypeName("variableSoftSphere");

    //- Construct from dictionary
    variableSoftSphere
    (
        const dictionary& dict,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~variableSoftSphere() = default;


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Return the collision cross section * relative velocity product
        virtual scalar sigmaTcR
        (
            const uspParcel& pP,
            const uspParcel& pQ
        ) const;

        //- Apply collision
        virtual void collide
        (
            uspParcel& pP,
            uspParcel& pQ,
            label& cellI
        );

        virtual const dictionary& propertiesDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "variableSoftSphere.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
