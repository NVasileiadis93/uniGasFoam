/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifndef gravitationalAccelerationController_H
#define gravitationalAccelerationController_H

#include "uspStateController.H"
#include "vector.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class gravitationalAccelerationController Declaration
\*---------------------------------------------------------------------------*/

class gravitationalAccelerationController
:
    public uspStateController
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;

        vector acceleration_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("gravitationalAccelerationController");

    //- Construct from components
    gravitationalAccelerationController
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~gravitationalAccelerationController() = default;

    // Member Functions

        virtual void initialConfiguration();
        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
