/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "LarsenBorgnakkeVariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphere, 0);
addToRunTimeSelectionTable
(
    BinaryCollisionModel,
    LarsenBorgnakkeVariableHardSphere,
    dictionary
);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::LarsenBorgnakkeVariableHardSphere
(
    const dictionary& dict,
    uspCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict("collisionCoeffs")),
    Tref_(coeffDict_.get<scalar>("Tref")),
    rotationalRelaxationCollisionNumber_
    (
        coeffDict_.get<scalar>("rotationalRelaxationCollisionNumber")
    ),
    electronicRelaxationCollisionNumber_
    (
        coeffDict_.get<scalar>("electronicRelaxationCollisionNumber")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LarsenBorgnakkeVariableHardSphere::active() const
{
    return true;
}


Foam::scalar Foam::LarsenBorgnakkeVariableHardSphere::sigmaTcR
(
    const uspParcel& pP,
    const uspParcel& pQ
) const
{
    const scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();

    const scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    const scalar mP = cloud_.constProps(typeIdP).mass();

    const scalar mQ = cloud_.constProps(typeIdQ).mass();

    const scalar mR = mP*mQ/(mP + mQ);

    // Calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    const scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


void Foam::LarsenBorgnakkeVariableHardSphere::collide
(
    uspParcel& pP,
    uspParcel& pQ,
    label& cellI
)
{
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    scalar& ERotP = pP.ERot();
    scalar& ERotQ = pQ.ERot();
    label& ELevelP = pP.ELevel();
    label& ELevelQ = pQ.ELevel();
    labelList& vibLevelP = pP.vibLevel();
    labelList& vibLevelQ = pQ.vibLevel();

    const vector positionP(pP.position());
    const vector positionQ(pQ.position());

    const scalar collisionSeparation =
        sqrt
        (
            sqr(positionP.x() - positionQ.x())
          + sqr(positionP.y() - positionQ.y())
        );

    cloud_.cellPropMeasurements().collisionSeparation()[cellI] +=
        collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;

    Random& rndGen(cloud_.rndGen());


    //   VIBRATIONAL ENERGY EXCHANGE - QUANTUM - KINETIC MODEL

    scalar preCollisionERotP = ERotP;
    scalar preCollisionERotQ = ERotQ;

    scalarList preCollisionEVibP(vibLevelP.size(), Zero);
    scalarList preCollisionEVibQ(vibLevelQ.size(), Zero);

    const label vibrationalDofP = cloud_.constProps(typeIdP).vibrationalDoF();
    const label vibrationalDofQ = cloud_.constProps(typeIdQ).vibrationalDoF();

    if (vibrationalDofP > 0)
    {
        forAll(vibLevelP, i)
        {
            preCollisionEVibP[i] =
                vibLevelP[i]
               *cloud_.constProps(typeIdP).thetaV()[i]
               *physicoChemical::k.value();
        }
    }

    if (vibrationalDofQ > 0)
    {
        forAll(vibLevelQ, i)
        {
            preCollisionEVibQ[i] =
                vibLevelQ[i]
               *cloud_.constProps(typeIdQ).thetaV()[i]
               *physicoChemical::k.value();
        }
    }

    const scalar preCollisionEEleP =
        cloud_.constProps(typeIdP).electronicEnergyList()[ELevelP];
    const scalar preCollisionEEleQ =
        cloud_.constProps(typeIdQ).electronicEnergyList()[ELevelQ];

    const label rotationalDofP = cloud_.constProps(typeIdP).rotationalDoF();
    const label rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDoF();


    label jMaxP = cloud_.constProps(typeIdP).nElectronicLevels();
    label jMaxQ = cloud_.constProps(typeIdQ).nElectronicLevels();

    const List<scalar>& EElistP =
        cloud_.constProps(typeIdP).electronicEnergyList();
    const List<scalar>& EElistQ =
        cloud_.constProps(typeIdQ).electronicEnergyList();

    const List<label>& gListP = cloud_.constProps(typeIdP).degeneracyList();
    const List<label>& gListQ = cloud_.constProps(typeIdQ).degeneracyList();

    const scalarList& thetaVP = cloud_.constProps(typeIdP).thetaV();
    const scalarList& thetaVQ = cloud_.constProps(typeIdQ).thetaV();

    const scalarList& thetaDP = cloud_.constProps(typeIdP).thetaD();
    const scalarList& thetaDQ = cloud_.constProps(typeIdQ).thetaD();

    const scalarList& ZrefP = cloud_.constProps(typeIdP).Zref();
    const scalarList& ZrefQ = cloud_.constProps(typeIdQ).Zref();

    const scalarList& refTempZvP = cloud_.constProps(typeIdP).TrefZv();
    const scalarList& refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar mP = cloud_.constProps(typeIdP).mass();
    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);
    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
    scalar cRsqr = magSqr(UP - UQ);

    scalar translationalEnergy = 0.5*mR*cRsqr;

    scalar ChiB = 2.5 - omegaPQ;

    scalar inverseRotationalCollisionNumber =
        1.0/rotationalRelaxationCollisionNumber_;
    scalar inverseElectronicCollisionNumber =
        1.0/electronicRelaxationCollisionNumber_;

//   Larsen Borgnakke rotational energy redistribution part.  Using the serial
//   application of the LB method, as per the INELRS subroutine in Bird's
//   USP0R.FOR

    if (inverseElectronicCollisionNumber > rndGen.sample01<scalar>())
    {
        // collision energy of particle P = relative translational energy +
        // pre-collision electronic energy

        scalar EcP = translationalEnergy + preCollisionEEleP;

        label postCollisionELevel =
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcP,
                jMaxP,
                omegaPQ,
                EElistP,
                gListP
            );

        ELevelP = postCollisionELevel;

        // Relative translational energy after electronic exchange
        translationalEnergy = EcP - EElistP[ELevelP];
    }

    if (vibrationalDofP > VSMALL)
    {
        forAll(vibLevelP, i)
        {
            // Collision energy of particle P = relative translational energy +
            // pre-collision vibrational energy
            scalar EcP = translationalEnergy + preCollisionEVibP[i];

            // Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxP = (EcP / (physicoChemical::k.value()*thetaVP[i]));

            if (iMaxP > SMALL)
            {
                vibLevelP[i] =
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        false,
                        vibLevelP[i],
                        iMaxP,
                        thetaVP[i],
                        thetaDP[i],
                        refTempZvP[i],
                        omegaPQ,
                        ZrefP[i],
                        EcP
                    );

                translationalEnergy =
                    EcP
                  - (
                        vibLevelP[i]
                       *cloud_.constProps(typeIdP).thetaV()[i]
                       *physicoChemical::k.value()
                    );
            }
        }
    }

    if (rotationalDofP > VSMALL)
    {
        if (inverseRotationalCollisionNumber > rndGen.sample01<scalar>())
        {
            scalar EcP = translationalEnergy + preCollisionERotP;

            scalar energyRatio =
                cloud_.postCollisionRotationalEnergy(rotationalDofP, ChiB);

            ERotP = energyRatio*EcP;

            translationalEnergy = EcP - ERotP;
        }
    }

    if (inverseElectronicCollisionNumber > rndGen.sample01<scalar>())
    {
        // Collision energy of particle Q = relative translational energy +
        // pre-collision electronic energy

        scalar EcQ = translationalEnergy + preCollisionEEleQ;

        const label postCollisionELevel =
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcQ,
                jMaxQ,
                omegaPQ,
                EElistQ,
                gListQ
            );

        ELevelQ = postCollisionELevel;

        // Relative translational energy after electronic exchange
        translationalEnergy = EcQ - EElistQ[ELevelQ];
    }

    if (vibrationalDofQ > VSMALL)
    {
        forAll(vibLevelQ, i)
        {
            // Collision energy of particle Q = relative translational energy +
            // pre-collision vibrational energy
            scalar EcQ = translationalEnergy + preCollisionEVibQ[i];

            // Maximum possible quantum level (equation 3, Bird 2010)
            label iMaxQ = (EcQ / (physicoChemical::k.value()*thetaVQ[i]));

            if (iMaxQ > SMALL)
            {
                vibLevelQ[i] =
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        false,
                        vibLevelQ[i],
                        iMaxQ,
                        thetaVQ[i],
                        thetaDQ[i],
                        refTempZvQ[i],
                        omegaPQ,
                        ZrefQ[i],
                        EcQ
                    );

                translationalEnergy =
                    EcQ
                  - (
                        vibLevelQ[i]
                       *cloud_.constProps(typeIdQ).thetaV()[i]
                       *physicoChemical::k.value()
                    );
            }
        }
    }

    if (rotationalDofQ > VSMALL)
    {
        if (inverseRotationalCollisionNumber > rndGen.sample01<scalar>())
        {
            scalar EcQ = translationalEnergy + preCollisionERotQ;

            scalar energyRatio =
                cloud_.postCollisionRotationalEnergy(rotationalDofQ, ChiB);

            ERotQ = energyRatio*EcQ;

            translationalEnergy = EcQ - ERotQ;
        }
    }

    // Rescale the translational energy
    const scalar cR = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part

    const scalar cosTheta = 2.0*rndGen.sample01<scalar>() - 1.0;

    const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    const scalar phi = 2.0*pi*rndGen.sample01<scalar>();

    const vector postCollisionRelU
    (
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        )
    );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
}


const Foam::dictionary&
Foam::LarsenBorgnakkeVariableHardSphere::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
