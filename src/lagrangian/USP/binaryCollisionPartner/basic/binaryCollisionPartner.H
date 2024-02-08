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
    Foam::binaryCollisionPartner

Description
    An abstract class for selecting collision partners in a USP simulation

SourceFiles
    binaryCollisionPartnerI.H
    binaryCollisionPartner.C

\*---------------------------------------------------------------------------*/

#ifndef binaryCollisionPartner_H
#define binaryCollisionPartner_H

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
                  Class binaryCollisionPartner Declaration
\*---------------------------------------------------------------------------*/

class binaryCollisionPartner
{
protected:

    // Protected data

        const fvMesh& mesh_;

        uspCloud& cloud_;

        Random& rndGen_;


public:

    //- Runtime type information
    TypeName("binaryCollisionPartner");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        binaryCollisionPartner,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uspCloud& cloud
        ),
        (dict, mesh, cloud)
    );

    //- Constructor
    binaryCollisionPartner
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Selectors
    static autoPtr<binaryCollisionPartner> New
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~binaryCollisionPartner() = default;


    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void collide() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "binaryCollisionPartnerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
