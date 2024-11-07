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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline List<autoPtr<uniGasPatchBoundary>>&
uniGasBoundaries::patchBoundaryModels()
{
    return patchBoundaryModels_;
}


inline const List<autoPtr<uniGasPatchBoundary>>&
uniGasBoundaries::patchBoundaryModels() const
{
    return patchBoundaryModels_;
}


inline const List<word>& uniGasBoundaries::patchBoundaryNames() const
{
    return patchBoundaryNames_;
}


inline const List<label>& uniGasBoundaries::patchBoundaryIds() const
{
    return patchBoundaryIds_;
}


inline const List<label>& uniGasBoundaries::patchToModelIds() const
{
    return patchToModelId_;
}


inline List<autoPtr<uniGasCyclicBoundary>>&
uniGasBoundaries::cyclicBoundaryModels()
{
    return cyclicBoundaryModels_;
}


inline const List<autoPtr<uniGasCyclicBoundary>>&
uniGasBoundaries::cyclicBoundaryModels() const
{
    return cyclicBoundaryModels_;
}


inline const List<word>& uniGasBoundaries::cyclicBoundaryNames() const
{
    return cyclicBoundaryNames_;
}


inline const List<label>& uniGasBoundaries::cyclicBoundaryIds() const
{
    return cyclicBoundaryIds_;
}


inline const List<label>& uniGasBoundaries::cyclicBoundaryToModelIds() const
{
    return cyclicBoundaryToModelId_;
}


inline List<autoPtr<uniGasGeneralBoundary>>&
uniGasBoundaries::generalBoundaryModels()
{
    return generalBoundaryModels_;
}


inline const List<autoPtr<uniGasGeneralBoundary>>&
uniGasBoundaries::generalBoundaryModels() const
{
    return generalBoundaryModels_;
}


inline const List<word>& uniGasBoundaries::generalBoundaryNames() const
{
    return generalBoundaryNames_;
}


inline const List<label>& uniGasBoundaries::generalBoundaryIds() const
{
    return generalBoundaryIds_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
