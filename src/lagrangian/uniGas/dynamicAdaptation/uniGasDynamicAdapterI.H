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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const label& uniGasDynamicAdapter::minSubCellLevels() const
{
    return minSubCellLevels_;
}


inline const label& uniGasDynamicAdapter::maxSubCellLevels() const
{
    return maxSubCellLevels_;
}


inline const bool& uniGasDynamicAdapter::cellWeightAdaptation() const
{
    return cellWeightAdaptation_;
}


inline const bool& uniGasDynamicAdapter::subCellAdaptation() const
{
    return subCellAdaptation_;
}


inline const scalar& uniGasDynamicAdapter::maxSubCellSizeMFPRatio() const
{
    return maxSubCellSizeMFPRatio_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
