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

template<class Type>
Foam::fvMesh& Foam::CloudWithModels<Type>::mesh() const
{
    return mesh_;
}


template<class Type>
const Foam::IOdictionary&
Foam::CloudWithModels<Type>::particleProperties() const
{
    return particleProperties_;
}


template<class Type>
const Foam::IOdictionary& Foam::CloudWithModels<Type>::outputProperties() const
{
    return outputProperties_;
}


template<class Type>
Foam::IOdictionary& Foam::CloudWithModels<Type>::outputProperties()
{
    return outputProperties_;
}


template<class Type>
const Foam::cloudSolution& Foam::CloudWithModels<Type>::solution() const
{
    return solution_;
}


template<class Type>
inline Foam::List<Foam::DynamicList<Type*>>&
Foam::CloudWithModels<Type>::cellOccupancy()
{
    if (!cellOccupancyPtr_)
    {
        buildCellOccupancy();
    }

    return *cellOccupancyPtr_;
}
