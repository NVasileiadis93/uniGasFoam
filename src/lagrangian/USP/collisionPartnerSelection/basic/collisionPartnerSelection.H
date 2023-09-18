/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::collisionPartnerSelection

Description
    An abstract class for selecting collision partners in a USP simulation

SourceFiles
    collisionPartnerSelectionI.H
    collisionPartnerSelection.C

\*---------------------------------------------------------------------------*/

#ifndef collisionPartnerSelection_H
#define collisionPartnerSelection_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class uspCloud;
class Random;

/*---------------------------------------------------------------------------*\
                  Class collisionPartnerSelection Declaration
\*---------------------------------------------------------------------------*/

class collisionPartnerSelection
{
protected:

    // Protected data

        const fvMesh& mesh_;

        uspCloud& cloud_;

        Random& rndGen_;


public:

    //- Runtime type information
    TypeName("collisionPartnerSelection");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        collisionPartnerSelection,
        dictionary,
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Constructor
    collisionPartnerSelection
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Selectors
    static autoPtr<collisionPartnerSelection> New
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~collisionPartnerSelection() = default;


    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void collide() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "collisionPartnerSelectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
