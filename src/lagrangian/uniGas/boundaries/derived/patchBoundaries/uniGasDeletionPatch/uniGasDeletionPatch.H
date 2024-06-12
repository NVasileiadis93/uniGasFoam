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
    Foam::uniGasDeletionPatch

Description

    Fully deletion boundary patch: uniGasParcels get deleted instantly upon
    collision with a face on the boundary (aka vacuum boundary condition).

    Use in the boundariesDict file.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch       outlet;
    }

    boundaryModel   uniGasDeletionPatch;

    uniGasDeletionPatchProperties
    {
    }
    \endverbatim

SourceFiles
    uniGasDeletionPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasDeletionPatch_H
#define uniGasDeletionPatch_H

#include "uniGasPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class uniGasDeletionPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasDeletionPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

        bool allSpecies_;

        labelList typeIds_;


    // Private Member Functions

public:

    //- Runtime type information
    TypeName("uniGasDeletionPatch");

    //- Construct from components
    uniGasDeletionPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasDeletionPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle
        (
            uniGasParcel& p,
            uniGasParcel::trackingData& td
        );

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
