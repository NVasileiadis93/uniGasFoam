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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

inline Foam::List<Foam::autoPtr<Foam::uniGasField>>&
Foam::uniGasFieldProperties::fields()
{
    return fields_;
}


inline const Foam::List<Foam::autoPtr<Foam::uniGasField>>&
Foam::uniGasFieldProperties::fields() const
{
    return fields_;
}


inline const Foam::List<Foam::word>&
Foam::uniGasFieldProperties::fieldNames() const
{
    return fieldNames_;
}


inline const Foam::List<Foam::label>&
Foam::uniGasFieldProperties::fieldIds() const
{
    return fieldIds_;
}


// ************************************************************************* //
