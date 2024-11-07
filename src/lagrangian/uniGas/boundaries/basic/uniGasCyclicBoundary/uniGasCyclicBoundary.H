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
    Foam::uniGasCyclicBoundary

Description
    Used for applying instantaneous actions on molecules colliding with a
    cyclic boundary.

SourceFiles
    uniGasCyclicBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasCyclicBoundary_H
#define uniGasCyclicBoundary_H

#include "autoPtr.H"
#include "uniGasBoundaryBase.H"
#include "uniGasParcel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class uniGasCyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class uniGasCyclicBoundary
:
    public uniGasBoundaryBase
{
protected:

    // Protected data

        //- Cyclic neighbour patch name
        word neighbPatchName_;

        //- Cyclic neighbour patch ID
        label neighbPatchId_;

        //- Faces on the uniGasCyclicBoundary
        labelList faces_;

        labelList coupledFacesA_;
        labelList coupledFacesB_;

        labelList cellsA_;
        labelList cellsB_;


    // Private Member Functions

        void getCoupledFaces(const polyMesh& mesh);


public:

    //- Runtime type information
    TypeName("uniGasCyclicBoundary");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uniGasCyclicBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    uniGasCyclicBoundary
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );


    //- Selector
    static autoPtr<uniGasCyclicBoundary> New
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~uniGasCyclicBoundary() = default;


    // Member Functions

        virtual void calculateProperties() = 0;

        //- Create an initial configuration
        virtual void initialConfiguration() = 0;

        //- Control the molecule
        virtual void controlMol
        (
            uniGasParcel& mol,
            uniGasParcel::trackingData& td
        ) = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;


        // Access

            //- Return faces on patch
            virtual const labelList& controlPatch() const;

            //- Return the adjacent cells
            virtual const labelList& controlZone() const;

            virtual const labelList& allFaces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
