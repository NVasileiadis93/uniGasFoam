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
    Foam::dsmcCollisionPartner

Description
    An abstract class for selecting collision partners for dsmc collision models

SourceFiles
    dsmcCollisionPartnerI.H
    dsmcCollisionPartner.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCollisionPartner_H
#define dsmcCollisionPartner_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class uniGasCloud;
class Random;

/*---------------------------------------------------------------------------*\
                  Class dsmcCollisionPartner Declaration
\*---------------------------------------------------------------------------*/

class dsmcCollisionPartner
{
protected:

    // Protected data

        const fvMesh& mesh_;

        uniGasCloud& cloud_;

        Random& rndGen_;


public:

    //- Runtime type information
    TypeName("dsmcCollisionPartner");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcCollisionPartner,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uniGasCloud& cloud
        ),
        (dict, mesh, cloud)
    );

    //- Constructor
    dsmcCollisionPartner
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& cloud
    );

    //- Selectors
    static autoPtr<dsmcCollisionPartner> New
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~dsmcCollisionPartner() = default;


    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void collide() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "dsmcCollisionPartnerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
