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

inline Foam::label Foam::uspCloud::binCollModel() const
{
    return binCollModel_;
}

inline Foam::label Foam::uspCloud::relCollModel() const
{
    return relCollModel_;
}

inline Foam::word Foam::uspCloud::relaxationCollisionModelName() const
{
    return relaxationCollisionModelName_;
}

inline Foam::boolVector Foam::uspCloud::solutionDimensions() const
{
    return solutionDimensions_;
}

inline const Foam::List<Foam::word>& Foam::uspCloud::typeIdList() const
{
    return typeIdList_;
}

inline Foam::scalar Foam::uspCloud::nParticle() const
{
    return nParticle_;
}

inline Foam::scalar Foam::uspCloud::collTref() const
{
    return collTref_;
}

inline bool Foam::uspCloud::adaptive() const
{
    return adaptive_;
}

inline bool Foam::uspCloud::cellWeighted() const
{
    return cellWeighted_;
}

inline Foam::label Foam::uspCloud::particlesPerSubcell() const
{
    return particlesPerSubcell_;
}

inline Foam::label Foam::uspCloud::minParticlesPerSubcell() const
{
    return minParticlesPerSubcell_;
}

inline Foam::scalar Foam::uspCloud::maxCellWeightRatio() const
{
    return maxCellWeightRatio_;
}

inline Foam::label Foam::uspCloud::maxSmoothingPasses() const
{
    return maxSmoothingPasses_;
}

inline Foam::scalar Foam::uspCloud::cellWF(const label& cell) const
{
    return cellWeightFactor_[cell];
}

inline Foam::label Foam::uspCloud::cellCollModel(const label& cell) const
{
    return label(cellCollisionModel_[cell]);
}

inline Foam::scalar Foam::uspCloud::radialExtent() const
{
    return radialExtent_;
}

inline Foam::scalar Foam::uspCloud::maxRWF() const
{
    return maxRWF_;
}

inline Foam::scalar Foam::uspCloud::axiRWF(const point& p) const
{
    if (axisymmetric_)
    {
        scalar radius = sqrt(p.y()*p.y() + p.z()*p.z());
        return 1.0 + (maxRWF_-1.0)*radius/radialExtent_;
    }

    return 1.0;
}

inline bool Foam::uspCloud::axisymmetric() const
{
    return axisymmetric_;
}

inline Foam::label Foam::uspCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}

inline Foam::volVectorField& Foam::uspCloud::subcellLevels()
{
    return subcellLevels_;
}

inline Foam::volScalarField& Foam::uspCloud::cellWeightFactor()
{
    return cellWeightFactor_;
}

inline Foam::volScalarField& Foam::uspCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

inline Foam::volScalarField& Foam::uspCloud::cellCollModel()
{
    return cellCollisionModel_;
}


inline Foam::scalarField& Foam::uspCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


inline const Foam::List<Foam::uspParcel::constantProperties>&
Foam::uspCloud::constProps() const
{
    return constProps_;
}


inline const Foam::uspParcel::constantProperties& Foam::uspCloud::constProps
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


inline Foam::Random& Foam::uspCloud::rndGen()
{
    return rndGen_;
}


inline const Foam::uspDynamicAdapter& Foam::uspCloud::dynamicAdapter() const
{
    return dynamicAdapter_;
}


inline Foam::uspDynamicAdapter& Foam::uspCloud::dynamicAdapter()
{
    return dynamicAdapter_;
}


inline const Foam::uspFieldProperties& Foam::uspCloud::fields() const
{
    return fields_;
}

inline Foam::uspFieldProperties& Foam::uspCloud::fields()
{
    return fields_;
}

inline const  Foam::uspBoundaries& Foam::uspCloud::boundaries() const
{
    return  boundaries_;
}


inline  Foam::uspBoundaries& Foam::uspCloud::boundaries()
{
    return  boundaries_;
}


inline const Foam::uspFaceTracker& Foam::uspCloud::tracker() const
{
    return trackingInfo_;
}


inline Foam::uspFaceTracker& Foam::uspCloud::tracker()
{
    return trackingInfo_;
}


inline const Foam::binaryCollisionModel& Foam::uspCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}


inline Foam::binaryCollisionModel& Foam::uspCloud::binaryCollision()
{
    return binaryCollisionModel_();
}


inline const Foam::relaxationCollisionModel& Foam::uspCloud::relaxationCollision() const
{
    return relaxationCollisionModel_;
}


inline Foam::relaxationCollisionModel& Foam::uspCloud::relaxationCollision()
{
    return relaxationCollisionModel_();
}


inline const Foam::uspHybridDecomposition& Foam::uspCloud::hybridDecomposition() const
{
    return hybridDecompositionModel_;
}


inline Foam::uspHybridDecomposition& Foam::uspCloud::hybridDecomposition()
{
    return hybridDecompositionModel_();
}


inline const Foam::uspReactions& Foam::uspCloud::reactions() const
{
    return reactions_;
}


inline Foam::uspReactions& Foam::uspCloud::reactions()
{
    return reactions_;
}


inline const Foam::boundaryMeasurements& Foam::uspCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}


inline Foam::boundaryMeasurements& Foam::uspCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}


inline const Foam::cellMeasurements& Foam::uspCloud::cellPropMeasurements() const
{
    return cellMeas_;
}


inline Foam::cellMeasurements& Foam::uspCloud::cellPropMeasurements()
{
    return cellMeas_;
}


inline Foam::scalarList Foam::uspCloud::infoMeasurements() const
{
    scalarList infoMeasurements(5, 0.0);

    forAllConstIter(uspCloud, *this, iter)
    {
        const uspParcel& p = iter();

        const uspParcel::constantProperties& cP = constProps(p.typeId());

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


inline Foam::scalar Foam::uspCloud::maxwellianAverageSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::tmp<Foam::scalarField> Foam::uspCloud::maxwellianAverageSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::scalar Foam::uspCloud::maxwellianRMSSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::tmp<Foam::scalarField> Foam::uspCloud::maxwellianRMSSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::scalar Foam::uspCloud::maxwellianMostProbableSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::tmp<Foam::scalarField> Foam::uspCloud::maxwellianMostProbableSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline void Foam::uspCloud::clear()
{
    return IDLList<uspParcel>::clear();
}

// ************************************************************************* //
