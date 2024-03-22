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

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const List<scalarField>& cellMeasurements::rhoNMean() const
{
    return rhoNMean_;
}


inline const List<scalarField>& cellMeasurements::rhoNInstantaneous() const
{
    return rhoNInstantaneous_;
}


inline const List<scalarField>& cellMeasurements::rhoNMeanXnParticle() const
{
    return rhoNMeanXnParticle_;
}


inline const List<scalarField>& cellMeasurements::rhoNMeanInt() const
{
    return rhoNMeanInt_;
}


inline const List<scalarField>& cellMeasurements::molsElec() const
{
    return molsElec_;
}


inline const List<scalarField>& cellMeasurements::rhoMMean() const
{
    return rhoMMean_;
}


inline const List<scalarField>& cellMeasurements::rhoMMeanXnParticle() const
{
    return rhoMMeanXnParticle_;
}


inline const List<scalarField>& cellMeasurements::linearKEMean() const
{
    return linearKEMean_;
}


inline const List<scalarField>& cellMeasurements::linearKEMeanXnParticle() const
{
    return linearKEMeanXnParticle_;
}


inline const List<scalarField>& cellMeasurements::rotationalEMean() const
{
    return rotationalEMean_;
}


inline const List<scalarField>& cellMeasurements::rotationalDofMean() const
{
    return rotationalDofMean_;
}


inline const List<scalarField>& cellMeasurements::muu() const
{
    return muu_;
}


inline const List<scalarField>& cellMeasurements::muv() const
{
    return muv_;
}


inline const List<scalarField>& cellMeasurements::muw() const
{
    return muw_;
}


inline const List<scalarField>& cellMeasurements::mvv() const
{
    return mvv_;
}


inline const List<scalarField>& cellMeasurements::mvw() const
{
    return mvw_;
}


inline const List<scalarField>& cellMeasurements::mww() const
{
    return mww_;
}


inline const List<scalarField>& cellMeasurements::mcc() const
{
    return mcc_;
}


inline const List<scalarField>& cellMeasurements::mccu() const
{
    return mccu_;
}


inline const List<scalarField>& cellMeasurements::mccv() const
{
    return mccv_;
}


inline const List<scalarField>& cellMeasurements::mccw() const
{
    return mccw_;
}


inline const List<scalarField>& cellMeasurements::eu() const
{
    return eu_;
}


inline const List<scalarField>& cellMeasurements::ev() const
{
    return ev_;
}


inline const List<scalarField>& cellMeasurements::ew() const
{
    return ew_;
}


inline const List<scalarField>& cellMeasurements::e() const
{
    return e_;
}


inline const List<vectorField>& cellMeasurements::momentumMean() const
{
    return momentumMean_;
}


inline const List<vectorField>& cellMeasurements::momentumMeanXnParticle() const
{
    return momentumMeanXnParticle_;
}


inline const List<List<scalarField>>& cellMeasurements::vibrationalETotal() const
{
    return vibrationalETotal_;
}


inline const List<scalarField>& cellMeasurements::electronicETotal() const
{
    return electronicETotal_;
}


inline const List<scalarField>& cellMeasurements::nParcels() const
{
    return nParcels_;
}


inline const List<scalarField>& cellMeasurements::nParcelsXnParticle() const
{
    return nParcelsXnParticle_;
}


inline const List<scalarField>& cellMeasurements::mccSpecies() const
{
    return mccSpecies_;
}


inline const List<scalarField>& cellMeasurements::nGroundElectronicLevel() const
{
    return nGroundElectronicLevel_;
}


inline const List<scalarField>& cellMeasurements::nFirstElectronicLevel() const
{
    return nFirstElectronicLevel_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
