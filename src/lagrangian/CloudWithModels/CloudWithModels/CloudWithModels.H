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
    Foam::CloudWithModels

Description
    Cloud with models class

SourceFiles
    CloudWithModelsI.H
    CloudWithModels.C

\*---------------------------------------------------------------------------*/

#ifndef CloudWithModels_H
#define CloudWithModels_H

#include "Cloud.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "cloudSolution.H"

namespace Foam
{
template<class Type>
class CloudWithModels
:
    public Cloud<Type>
{
public:

    // Public typedefs

        //- Type of particle the cloud was instantiated for
        typedef Type parcelType;

        typedef CloudWithModels<Type> baseCloudType;


protected:

    // Protected data

        //- Reference to the fvMesh
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Dictionary of output properties
        IOdictionary outputProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Cell occupancy information for each parcel, (demand driven)
        autoPtr<List<DynamicList<Type*>>> cellOccupancyPtr_;


    // Protected Member Functions


public:

    // Constructors

        //- Construct from mesh and a list of particles
        CloudWithModels
        (
            const fvMesh& mesh,
            const word& cloudName,
            const IDLList<Type>& particles
        );

        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        CloudWithModels
        (
            const fvMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        //- Return a reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;

        //- Return output properties dictionary
        inline const IOdictionary& outputProperties() const;

        //- Return non-const access to the output properties dictionary
        inline IOdictionary& outputProperties();

        //- Return const access to the solution properties
        inline const cloudSolution& solution() const;


    // Cell occupancy

        //- Return the cell occupancy information for each parcel, non-const
        //- access, the caller is responsible for updating it for its own
        //- purposes if particles are removed or created.
        inline List<DynamicList<parcelType*>>& cellOccupancy();

        //- Build the cellOccupancy
        void buildCellOccupancy();

        //- Update/build the cellOccupancy if it has already been used
        void updateCellOccupancy();

        //- Insert parcel
        void insertParcelInCellOccupancy(parcelType* p);

        //- Remove parcel
        void removeParcelFromCellOccupancy
        (
            const label parceli,
            const label celli
        );
};

} // End namespace Foam

#ifdef NoRepository
    #include "CloudWithModels.C"
#endif

#include "CloudWithModelsI.H"

#endif
