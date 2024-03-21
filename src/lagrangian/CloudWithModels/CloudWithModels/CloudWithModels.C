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

\*---------------------------------------------------------------------------*/

#include "CloudWithModels.H"

template<class Type>
Foam::CloudWithModels<Type>::CloudWithModels
(
    fvMesh& mesh,
    const word& cloudName,
    const IDLList<Type>& particles
)
:
    Cloud<Type>(mesh, cloudName, particles),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    cellOccupancyPtr_(nullptr)
{}


template<class Type>
Foam::CloudWithModels<Type>::CloudWithModels
(
    fvMesh& mesh,
    const word& cloudName,
    const bool checkClass
)
:
    Cloud<Type>(mesh, cloudName, checkClass),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    outputProperties_
    (
        IOobject
        (
            cloudName + "OutputProperties",
            mesh.time().timeName(),
            "uniform"/cloud::prefix/cloudName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    solution_(mesh_, particleProperties_.subDict("solution")),
    cellOccupancyPtr_(nullptr)
{}


template<class Type>
void Foam::CloudWithModels<Type>::buildCellOccupancy()
{
    if (!cellOccupancyPtr_)
    {
        cellOccupancyPtr_.reset
        (
            new List<DynamicList<parcelType*>>(mesh_.nCells())
        );
    }
    else if (cellOccupancyPtr_().size() != mesh_.nCells())
    {
        // If the size of the mesh has changed, reset the
        // cellOccupancy size

        cellOccupancyPtr_().setSize(mesh_.nCells());
    }

    auto& cellOccupancy = cellOccupancyPtr_();

    for (auto& list : cellOccupancy)
    {
        list.clear();
    }

    for (parcelType& p : *this)
    {
        cellOccupancy[p.cell()].append(&p);
    }
}


template<class Type>
void Foam::CloudWithModels<Type>::updateCellOccupancy()
{
    // Only build the cellOccupancy if the pointer is set, i.e. it has
    // been requested before.

    if (cellOccupancyPtr_)
    {
        buildCellOccupancy();
    }
}


template<class Type>
void Foam::CloudWithModels<Type>::insertParcelInCellOccupancy(parcelType* p)
{
    if (cellOccupancyPtr_)
    {
        auto& cellOccupancy = cellOccupancyPtr_();
        cellOccupancy[p->cell()].append(p);
    }
}


template<class Type>
void Foam::CloudWithModels<Type>::removeParcelFromCellOccupancy
(
    const label parceli,
    const label celli
)
{
    if (cellOccupancyPtr_)
    {
        auto& cellOccupancy = cellOccupancyPtr_();

        DynamicList<parcelType*> occupancy(cellOccupancy[celli].size());

        forAll(cellOccupancy[celli], pi)
        {
            if (pi != parceli)
            {
                occupancy.append(cellOccupancy[celli][pi]);
            }
        }

        cellOccupancy[celli].clear();
        cellOccupancy[celli].transfer(occupancy);
    }
}
