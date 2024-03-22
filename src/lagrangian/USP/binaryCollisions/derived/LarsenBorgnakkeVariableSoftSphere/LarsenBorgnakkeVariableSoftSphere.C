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

#include "LarsenBorgnakkeVariableSoftSphere.H"
#include "constants.H"
#include "uspCloud.H"
#include "uspParcel.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(LarsenBorgnakkeVariableSoftSphere, 0);
addToRunTimeSelectionTable
(
    binaryCollisionModel,
    LarsenBorgnakkeVariableSoftSphere,
    dictionary
);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableSoftSphere::LarsenBorgnakkeVariableSoftSphere
(
    const dictionary& dict,
    uspCloud& cloud
)
:
    binaryCollisionModel(dict, cloud),
    propertiesDict_(dict.subDict("collisionProperties")),
    Tref_(propertiesDict_.get<scalar>("Tref")),
    rotationalRelaxationCollisionNumber_
    (
        propertiesDict_.get<scalar>("rotationalRelaxationCollisionNumber")
    ),
    electronicRelaxationCollisionNumber_
    (
        propertiesDict_.get<scalar>("electronicRelaxationCollisionNumber")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LarsenBorgnakkeVariableSoftSphere::active() const
{
    return true;
}


Foam::scalar Foam::LarsenBorgnakkeVariableSoftSphere::sigmaTcR
(
    const uspParcel& pP,
    const uspParcel& pQ
) const
{

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


void Foam::LarsenBorgnakkeVariableSoftSphere::collide
(
    uspParcel& pP,
    uspParcel& pQ,
    label& cellI
)
{
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    scalar& ERotP = pP.ERot();
    scalar& ERotQ = pQ.ERot();
    label& ELevelP = pP.ELevel();
    label& ELevelQ = pQ.ELevel();
    labelList& vibLevelP = pP.vibLevel();
    labelList& vibLevelQ = pQ.vibLevel();

    const label typeIdP = pP.typeId();
    const uspParcel::constantProperties& cpP = cloud_.constProps(typeIdP);
    const label typeIdQ = pQ.typeId();
    const uspParcel::constantProperties& cpQ = cloud_.constProps(typeIdQ);

    scalar alphaPQ = 0.5*(cpP.alpha() + cpQ.alpha());

    Random& rndGen(cloud_.rndGen());

    // Larsen Borgnakke rotational energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // USP0R.FOR

    scalar preCollisionERotP = ERotP;
    scalar preCollisionERotQ = ERotQ;

    scalarList preCollisionEVibP(vibLevelP.size(), 0.0);
    scalarList preCollisionEVibQ(vibLevelQ.size(), 0.0);

    label vibrationalDofP = cpP.vibrationalDoF();
    label vibrationalDofQ = cpQ.vibrationalDoF();

    if (vibrationalDofP > 0)
    {
        forAll(vibLevelP, i)
        {
            preCollisionEVibP[i] =
                vibLevelP[i]
               *cpP.thetaV()[i]
               *physicoChemical::k.value();
        }
    }

    if (vibrationalDofQ > 0)
    {
        forAll(vibLevelQ, i)
        {
            preCollisionEVibQ[i] =
                vibLevelQ[i]
               *cpQ.thetaV()[i]
               *physicoChemical::k.value();
        }
    }

    scalar preCollisionEEleP = cpP.electronicEnergyList()[ELevelP];
    scalar preCollisionEEleQ = cpQ.electronicEnergyList()[ELevelQ];

    label rotationalDofP = cpP.rotationalDoF();
    label rotationalDofQ = cpQ.rotationalDoF();

    label jMaxP = cpP.nElectronicLevels();
    label jMaxQ = cpQ.nElectronicLevels();

    const List<scalar>& EElistP = cpP.electronicEnergyList();
    const List<scalar>& EElistQ = cpQ.electronicEnergyList();

    const List<label>& gListP = cpP.degeneracyList();
    const List<label>& gListQ = cpQ.degeneracyList();

    const scalarList& thetaVP = cpP.thetaV();
    const scalarList& thetaVQ = cpQ.thetaV();

    const scalarList& thetaDP = cpP.thetaD();
    const scalarList& thetaDQ = cpQ.thetaD();

    const scalarList& ZrefP = cpP.Zref();
    const scalarList& ZrefQ = cpQ.Zref();

    const scalarList& refTempZvP = cpP.TrefZv();
    const scalarList& refTempZvQ = cpQ.TrefZv();

    scalar omegaPQ = 0.5*(cpP.omega() + cpQ.omega());

    scalar mP = cpP.mass();
    scalar mQ = cpQ.mass();

    scalar mR = mP*mQ/(mP + mQ);
    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
    scalar cRsqr = magSqr(UP - UQ);

    scalar translationalEnergy = 0.5*mR*cRsqr;

    scalar ChiB = 2.5 - omegaPQ;

    scalar inverseRotationalCollisionNumber =
        1.0/rotationalRelaxationCollisionNumber_;
    scalar inverseElectronicCollisionNumber =
        1.0/electronicRelaxationCollisionNumber_;

    // Larsen Borgnakke rotational energy redistribution part.  Using the
    // serial application of the LB method, as per the INELRS subroutine in
    // Bird's USP0R.FOR

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

        // relative translational energy after electronic exchange
        translationalEnergy = EcP - EElistP[ELevelP];
    }

    if (vibrationalDofP > VSMALL)
    {
        forAll(vibLevelP, i)
        {
            // collision energy of particle P = relative translational energy +
            // pre-collision vibrational energy
            scalar EcP = translationalEnergy + preCollisionEVibP[i];

            // maximum possible quantum level (equation 3, Bird 2010)
            label iMaxP = (EcP / (physicoChemical::k.value()*thetaVP[i]));

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
                       *physicoChemical::k.value()
                    );
            }
        }
    }

    if (rotationalDofP > 0)
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
        // collision energy of particle Q = relative translational energy +
        // pre-collision electronic energy

        scalar EcQ = translationalEnergy + preCollisionEEleQ;

        label postCollisionELevel =
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcQ,
                jMaxQ,
                omegaPQ,
                EElistQ,
                gListQ
            );

        ELevelQ = postCollisionELevel;

        // relative translational energy after electronic exchange
        translationalEnergy = EcQ - EElistQ[ELevelQ];
    }

    if (vibrationalDofQ > VSMALL)
    {
        forAll(vibLevelQ, i)
        {
            // collision energy of particle Q = relative translational energy +
            // pre-collision vibrational energy
            scalar EcQ = translationalEnergy + preCollisionEVibQ[i];

            // - maximum possible quantum level (equation 3, Bird 2010)
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
                       *cpQ.thetaV()[i]
                       *physicoChemical::k.value()
                    );
            }
        }
    }

    if (rotationalDofQ > 0)
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
    scalar A = sqrt(2.0*translationalEnergy/mR);

    scalar cR = mag(UP - UQ);

    vector cRComponents = (UP - UQ)*(A/cR);

    cR = A;

    scalar cosTheta = (2.0*pow(rndGen.sample01<scalar>(), (1.0/alphaPQ))) - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.sample01<scalar>();

    scalar D = sqrt(pow(cRComponents.y(), 2.0) + pow(cRComponents.z(), 2.0));

    vector postCollisionRelU =
        vector
        (
            (cosTheta*cRComponents.x()) + sinTheta*sin(phi)*D,
            (cosTheta*cRComponents.y()) + sinTheta
               *(
                   cR*cRComponents.z()*cos(phi)
                 - cRComponents.x()*cRComponents.y()*sin(phi)
                )/D,
            (cosTheta*cRComponents.z())
              - sinTheta
               *(
                    cR*cRComponents.y()*cos(phi)
                  + cRComponents.x()*cRComponents.z()*sin(phi)
                )/D
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);

}


const Foam::dictionary&
Foam::LarsenBorgnakkeVariableSoftSphere::propertiesDict() const
{
    return propertiesDict_;
}


// ************************************************************************* //
