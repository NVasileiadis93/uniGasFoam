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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::uniGasCloud::binCollModel() const
{
    return binCollModel_;
}

inline Foam::label Foam::uniGasCloud::relCollModel() const
{
    return relCollModel_;
}

inline Foam::word Foam::uniGasCloud::relaxationCollisionModelName() const
{
    return relaxationCollisionModelName_;
}

inline Foam::boolVector Foam::uniGasCloud::solutionDimensions() const
{
    return solutionDimensions_;
}

inline const Foam::List<Foam::word>& Foam::uniGasCloud::typeIdList() const
{
    return typeIdList_;
}

inline Foam::scalar Foam::uniGasCloud::nParticle() const
{
    return nParticle_;
}

inline Foam::scalar Foam::uniGasCloud::collTref() const
{
    return collTref_;
}

inline bool Foam::uniGasCloud::adaptive() const
{
    return adaptive_;
}

inline bool Foam::uniGasCloud::cellWeighted() const
{
    return cellWeighted_;
}

inline Foam::label Foam::uniGasCloud::particlesPerSubCell() const
{
    return particlesPerSubCell_;
}

inline Foam::label Foam::uniGasCloud::minParticlesPerSubCell() const
{
    return minParticlesPerSubCell_;
}

inline Foam::scalar Foam::uniGasCloud::maxCellWeightRatio() const
{
    return maxCellWeightRatio_;
}

inline Foam::label Foam::uniGasCloud::maxSmoothingPasses() const
{
    return maxSmoothingPasses_;
}

inline Foam::scalar Foam::uniGasCloud::cellWF(const label& cell) const
{
    return cellWeightFactor_[cell];
}

inline Foam::label Foam::uniGasCloud::cellCollModel(const label& cell) const
{
    return label(cellCollisionModel_[cell]);
}

inline Foam::scalar Foam::uniGasCloud::radialExtent() const
{
    return radialExtent_;
}

inline Foam::scalar Foam::uniGasCloud::maxRWF() const
{
    return maxRWF_;
}

inline Foam::scalar Foam::uniGasCloud::axiRWF(const point& p) const
{
    if (axisymmetric_)
    {
        scalar radius = sqrt(p.y()*p.y() + p.z()*p.z());
        return 1.0 + (maxRWF_-1.0)*radius/radialExtent_;
    }

    return 1.0;
}

inline bool Foam::uniGasCloud::axisymmetric() const
{
    return axisymmetric_;
}

inline Foam::label Foam::uniGasCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}

inline Foam::volVectorField& Foam::uniGasCloud::subCellLevels()
{
    return subCellLevels_;
}

inline Foam::volScalarField& Foam::uniGasCloud::cellWeightFactor()
{
    return cellWeightFactor_;
}

inline Foam::volScalarField& Foam::uniGasCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

inline Foam::volScalarField& Foam::uniGasCloud::cellCollModel()
{
    return cellCollisionModel_;
}


inline Foam::scalarField& Foam::uniGasCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


inline const Foam::List<Foam::uniGasParcel::constantProperties>&
Foam::uniGasCloud::constProps() const
{
    return constProps_;
}


inline const Foam::uniGasParcel::constantProperties& Foam::uniGasCloud::constProps
(
    const label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


inline Foam::Random& Foam::uniGasCloud::rndGen()
{
    return rndGen_;
}


inline const Foam::uniGasDynamicAdapter& Foam::uniGasCloud::dynamicAdapter() const
{
    return dynamicAdapter_;
}


inline Foam::uniGasDynamicAdapter& Foam::uniGasCloud::dynamicAdapter()
{
    return dynamicAdapter_;
}


inline const Foam::uniGasFieldProperties& Foam::uniGasCloud::fields() const
{
    return fields_;
}

inline Foam::uniGasFieldProperties& Foam::uniGasCloud::fields()
{
    return fields_;
}

inline const  Foam::uniGasBoundaries& Foam::uniGasCloud::boundaries() const
{
    return  boundaries_;
}


inline  Foam::uniGasBoundaries& Foam::uniGasCloud::boundaries()
{
    return  boundaries_;
}


inline const Foam::uniGasFaceTracker& Foam::uniGasCloud::tracker() const
{
    return trackingInfo_;
}


inline Foam::uniGasFaceTracker& Foam::uniGasCloud::tracker()
{
    return trackingInfo_;
}


inline const Foam::binaryCollisionModel& Foam::uniGasCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}


inline Foam::binaryCollisionModel& Foam::uniGasCloud::binaryCollision()
{
    return binaryCollisionModel_();
}


inline const Foam::relaxationCollisionModel& Foam::uniGasCloud::relaxationCollision() const
{
    return relaxationCollisionModel_;
}


inline Foam::relaxationCollisionModel& Foam::uniGasCloud::relaxationCollision()
{
    return relaxationCollisionModel_();
}


inline const Foam::uniGasHybridDecomposition& Foam::uniGasCloud::hybridDecomposition() const
{
    return hybridDecompositionModel_;
}


inline Foam::uniGasHybridDecomposition& Foam::uniGasCloud::hybridDecomposition()
{
    return hybridDecompositionModel_();
}


inline const Foam::uniGasReactions& Foam::uniGasCloud::reactions() const
{
    return reactions_;
}


inline Foam::uniGasReactions& Foam::uniGasCloud::reactions()
{
    return reactions_;
}


inline const Foam::boundaryMeasurements& Foam::uniGasCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}


inline Foam::boundaryMeasurements& Foam::uniGasCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}


inline const Foam::cellMeasurements& Foam::uniGasCloud::cellPropMeasurements() const
{
    return cellMeas_;
}


inline Foam::cellMeasurements& Foam::uniGasCloud::cellPropMeasurements()
{
    return cellMeas_;
}


inline Foam::scalarList Foam::uniGasCloud::infoMeasurements() const
{
    scalarList infoMeasurements(5, 0.0);

    forAllConstIter(uniGasCloud, *this, iter)
    {
        const uniGasParcel& p = iter();

        const uniGasParcel::constantProperties& cP = constProps(p.typeId());

        // massInSystem
        infoMeasurements[0] += cP.mass();

        // linearKineticEnergy
        infoMeasurements[1] += 0.5*cP.mass()*(p.U() & p.U());

        // rotationalEnergy
        infoMeasurements[2] += p.ERot();
        if (cP.vibrationalDoF() > 0)
        {
            forAll(p.vibLevel(), i)
            {
                // vibrationalEnergy
                infoMeasurements[3] +=
                    physicoChemical::k.value()
                   *p.vibLevel()[i]*cP.thetaV()[i];
            }
        }

        // electronicEnergy
        infoMeasurements[4] += cP.electronicEnergyList()[p.ELevel()];
    }

    infoMeasurements[0] *= nParticle_;
    infoMeasurements[1] *= nParticle_;
    infoMeasurements[2] *= nParticle_;
    infoMeasurements[3] *= nParticle_;
    infoMeasurements[4] *= nParticle_;

    return infoMeasurements;
}


inline Foam::scalar Foam::uniGasCloud::maxwellianAverageSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::tmp<Foam::scalarField> Foam::uniGasCloud::maxwellianAverageSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::scalar Foam::uniGasCloud::maxwellianRMSSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::tmp<Foam::scalarField> Foam::uniGasCloud::maxwellianRMSSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::scalar Foam::uniGasCloud::maxwellianMostProbableSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::tmp<Foam::scalarField> Foam::uniGasCloud::maxwellianMostProbableSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline void Foam::uniGasCloud::clear()
{
    return IDLList<uniGasParcel>::clear();
}

// ************************************************************************* //
