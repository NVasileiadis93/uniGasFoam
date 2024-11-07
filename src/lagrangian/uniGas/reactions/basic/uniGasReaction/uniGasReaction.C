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

#include "uniGasReaction.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasReaction, 0);
defineRunTimeSelectionTable(uniGasReaction, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::uniGasReaction::dissociatePQ
(
    const scalar& heatOfReaction,
    const labelPair& productIDs,
    uniGasParcel& p,
    uniGasParcel& q
) const
{
    const scalar ChiB = 2.5 - this->omega(p, q);

    const scalar heatOfReactionDissociationJoules =
        heatOfReaction*physicoChemical::k.value();

    scalar translationalEnergy = this->translationalEnergy(p, q);

    translationalEnergy += heatOfReactionDissociationJoules + this->EVib(p);

    translationalEnergy += this->EEle(q);

    const label ELevelQ =
        cloud_.postCollisionElectronicEnergyLevel
        (
            translationalEnergy,
            this->jMax(q) + 1,
            this->omega(p, q),
            this->EEList(q),
            this->gList(q)
        );

    translationalEnergy -= this->EEList(q)[ELevelQ];

    label vibLevelQ = 0;
    scalar ERotQ = 0.0;

    if (this->rotationalDof(q) > VSMALL)
    {
        translationalEnergy += this->EVib(q);

        label iMax =
            translationalEnergy
           /(physicoChemical::k.value()*this->thetaV(q));

        vibLevelQ =
            cloud_.postCollisionVibrationalEnergyLevel
            (
                true,
                q.vibLevel()[0],
                iMax,
                this->thetaV(q),
                this->thetaD(q),
                this->refTempZv(q),
                this->omega(p, q),
                this->Zref(q),
                translationalEnergy
            );

        translationalEnergy -=
            vibLevelQ*this->thetaV(q)*physicoChemical::k.value();

        translationalEnergy += this->ERot(q);

        ERotQ =
            translationalEnergy
           *cloud_.postCollisionRotationalEnergy(this->rotationalDof(q), ChiB);

        translationalEnergy -= ERotQ;
    }

    const scalar relVel = sqrt(2.0*translationalEnergy/this->mR(p, q));

    // Variable Hard Sphere collision part

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelU =
        relVel
        *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    const vector UP =
        this->Ucm(p, q)
      + postCollisionRelU*this->m(q)/(this->m(p) + this->m(q));
    const vector UQ =
        this->Ucm(p, q)
      - postCollisionRelU*this->m(p)/(this->m(p) + this->m(q));

    // Mass of Product one and two
    const scalar mP1 = cloud_.constProps(productIDs[0]).mass();
    const scalar mP2 = cloud_.constProps(productIDs[1]).mass();

    const scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERot(p) + this->EEle(p);

    const scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    const scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    const scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    const scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelU2 =
        cRatoms
       *vector
        (
            cosTheta2,
            sinTheta2*cos(phi2),
            sinTheta2*sin(phi2)
        );

    const vector UP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
    const vector UP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

    q.U() = UQ;
    q.ERot() = ERotQ;
    if (this->rotationalDof(q) > VSMALL)
    {
        q.vibLevel()[0] = vibLevelQ;
    }

    q.ELevel() = ELevelQ;

    // Molecule P will dissociate
    const vector position(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    const scalar CWF = p.CWF();
    const scalar RWF = p.RWF();

    p.position() = position;
    p.U() = UP1;
    p.CWF() = CWF;
    p.RWF() = RWF;
    p.ERot() = 0.0;
    p.ELevel() = 0;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = productIDs[0];
    p.newParcel() = 0;
    p.vibLevel() = labelList();

    // Insert new product
    cloud_.addNewParcel
    (
        position,
        UP2,
        CWF,
        RWF,
        0.0,
        0,
        cell,
        productIDs[1],
        0,
        labelList()
    );
}


void Foam::uniGasReaction::ionisePQ
(
    const scalar& heatOfReaction,
    const labelPair& productIDs,
    uniGasParcel& p,
    uniGasParcel& q
) const
{
    const scalar ChiB = 2.5 - this->omega(p, q);

    const scalar heatOfReactionIonisationJoules =
        heatOfReaction*physicoChemical::k.value();

    scalar translationalEnergy = this->translationalEnergy(p, q);

    translationalEnergy += heatOfReactionIonisationJoules + this->EEle(p);

    translationalEnergy += this->EEle(q);

    label ELevelQ =
        cloud_.postCollisionElectronicEnergyLevel
        (
            translationalEnergy,
            this->jMax(q) + 1,
            this->omega(p, q),
            this->EEList(q),
            this->gList(q)
        );

    translationalEnergy -= this->EEList(q)[ELevelQ];

    label vibLevelQ = 0;
    scalar ERotQ = 0.0;

    if (this->rotationalDof(q) > VSMALL)
    {
        translationalEnergy += this->EVib(q);

        label iMax =
            translationalEnergy/(physicoChemical::k.value()*this->thetaV(q));

        vibLevelQ =
            cloud_.postCollisionVibrationalEnergyLevel
            (
                true,
                q.vibLevel()[0],
                iMax,
                this->thetaV(q),
                this->thetaD(q),
                this->refTempZv(q),
                this->omega(p, q),
                this->Zref(q),
                translationalEnergy
            );

        translationalEnergy -=
            vibLevelQ*this->thetaV(q)*physicoChemical::k.value();

        translationalEnergy += this->ERot(q);

        ERotQ =
            translationalEnergy
           *cloud_.postCollisionRotationalEnergy(this->rotationalDof(q), ChiB);

        translationalEnergy -= ERotQ;
    }

    const scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/this->mR(p, q));

    // Variable Hard Sphere collision part

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelU =
        relVelNonDissoMol
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    const vector UP =
        this->Ucm(p, q)
      + postCollisionRelU*this->m(q)/(this->m(p) + this->m(q));
    const vector UQ =
        this->Ucm(p, q)
      - postCollisionRelU*this->m(p)/(this->m(p) + this->m(q));

    // Mass of Product one and two
    const scalar mP1 = cloud_.constProps(productIDs[0]).mass();
    const scalar mP2 = cloud_.constProps(productIDs[1]).mass();

    const scalar mRatoms = mP1*mP2/(mP1 + mP2);

    translationalEnergy = this->ERot(p) + this->EVib(p);

    const scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

    // Variable Hard Sphere collision part
    const scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    const scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);

    const scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelU2 =
        cRatoms
       *vector
        (
            cosTheta2,
            sinTheta2*cos(phi2),
            sinTheta2*sin(phi2)
        );

    const vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
    const vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

    q.U() = UQ;
    q.ERot() = ERotQ;
    if (this->rotationalDof(q) > VSMALL)
    {
        q.vibLevel()[0] = vibLevelQ;
    }
    q.ELevel() = ELevelQ;

    const vector position(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );

    const scalar CWF = p.CWF();
    const scalar RWF = p.RWF();
    labelList vibLevel(0, 0);
    labelList electronVibLevel(1, 0);
    electronVibLevel[0] = 0;

    if (this->rotationalDof(p) > VSMALL)
    {
        vibLevel.resize(1);
        vibLevel[0] = 0;
    }

    p.position() = position;
    p.U() = uP1;
    p.CWF() = CWF;
    p.RWF() = RWF;
    p.ERot() = 0.0;
    p.ELevel() = 0;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = productIDs[0];
    p.newParcel() = 0;
    p.vibLevel() = vibLevel;

    // Insert new product (electron)
    cloud_.addNewParcel
    (
        position,
        uP2,
        CWF,
        RWF,
        0.0,
        0,
        cell,
        productIDs[1],
        0,
        electronVibLevel
    );
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::uniGasReaction::EVib
(
    const uniGasParcel& p
) const
{
    const label id = p.typeId();

    if (cloud_.constProps(id).vibrationalDoF() > VSMALL)
    {
        return
            p.vibLevel()[0]
           *cloud_.constProps(id).thetaV()[0]
           *physicoChemical::k.value();
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar Foam::uniGasReaction::ERot
(
    const uniGasParcel& p
) const
{
    return p.ERot();
}


Foam::scalar Foam::uniGasReaction::EEle
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).electronicEnergyList()[p.ELevel()];
}


Foam::scalar Foam::uniGasReaction::m
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).mass();
}


Foam::scalar Foam::uniGasReaction::mR
(
    const uniGasParcel& p,
    const uniGasParcel& q
) const
{
    return m(p)*m(q)/(m(p) + m(q));
}


Foam::scalar Foam::uniGasReaction::translationalEnergy
(
    const uniGasParcel& p,
    const uniGasParcel& q
) const
{
    return 0.5*mR(p, q)*magSqr(p.U() - q.U());
}


Foam::scalar Foam::uniGasReaction::omega
(
    const uniGasParcel& p,
    const uniGasParcel& q
) const
{
    return
        0.5
       *(
            cloud_.constProps(p.typeId()).omega()
          + cloud_.constProps(q.typeId()).omega()
        );
}



Foam::scalar Foam::uniGasReaction::thetaD
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).thetaD()[0];
}


Foam::scalar Foam::uniGasReaction::thetaV
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).thetaV()[0];
}


Foam::scalar Foam::uniGasReaction::Zref
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).Zref()[0];
}


Foam::scalar Foam::uniGasReaction::refTempZv
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).TrefZv()[0];
}


Foam::label Foam::uniGasReaction::charDissLevel
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).charDissQuantumLevel()[0];
}


Foam::label Foam::uniGasReaction::jMax
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).nElectronicLevels() - 1;
}


Foam::label Foam::uniGasReaction::rotationalDof
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).rotationalDoF();
}


Foam::vector Foam::uniGasReaction::Ucm
(
    const uniGasParcel& p,
    const uniGasParcel& q
) const
{
    scalar mP = m(p);
    scalar mQ = m(q);

    return (mP*p.U() + mQ*q.U())/(mP + mQ);
}


Foam::scalarList Foam::uniGasReaction::EEList
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).electronicEnergyList();
}



Foam::labelList Foam::uniGasReaction::gList
(
    const uniGasParcel& p
) const
{
    return cloud_.constProps(p.typeId()).degeneracyList();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasReaction::uniGasReaction
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    nTotReactions_(0),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactionName_(propsDict_.get<word>("reaction")),
    reactants_(),
    reactantIds_(),
    rDof1_(),
    rDof2_(),
    vDof1_(),
    vDof2_(),
    charge1_(),
    charge2_(),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uniGasReaction> Foam::uniGasReaction::New
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
{
    word modelType(dict.get<word>("reactionModel"));

    Info<< "Selecting the reaction model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "reactionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uniGasReaction>(cstrIter()(t, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasReaction::setCommonReactionProperties()
{
    // Reading in reactants

    reactants_ = propsDict_.get<wordPair>("reactants");

    reactantIds_ = labelPair(-1, -1);

    allowSplitting_ = propsDict_.get<bool>("allowSplitting");

    writeRatesToTerminal_ = propsDict_.get<bool>("writeRatesToTerminal");

    forAll(reactants_, i)
    {
        const word& reactantName = reactants_[i];

        reactantIds_[i] = cloud_.typeIdList().find(reactantName);

        // check that reactants belong to the typeIdList
        if (reactantIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find reactant " << reactantName << nl
                << "Available reactants:" << cloud_.typeIdList()
                << exit(FatalError);
        }
    }

    rDof1_ = cloud_.constProps(reactantIds_[0]).rotationalDoF();
    rDof2_ = cloud_.constProps(reactantIds_[1]).rotationalDoF();

    vDof1_ = cloud_.constProps(reactantIds_[0]).vibrationalDoF();
    vDof2_ = cloud_.constProps(reactantIds_[1]).vibrationalDoF();

    charge1_ = cloud_.constProps(reactantIds_[0]).charge();
    charge2_ = cloud_.constProps(reactantIds_[1]).charge();
}


void Foam::uniGasReaction::associativeIonisation
(
    const scalar& heatOfReactionIntermediateIonisation,
    const scalar& heatOfReactionRecombination,
    const labelPair& assIonProductIds,
    uniGasParcel& p,
    uniGasParcel& q
)
{
    scalar heatOfReactionIonisationJoules =
        heatOfReactionIntermediateIonisation
       *physicoChemical::k.value();
    scalar heatOfReactionRecombinationJoules =
        heatOfReactionRecombination
       *physicoChemical::k.value();

    scalar translationalEnergy = this->translationalEnergy(p, q);

    translationalEnergy +=
        heatOfReactionRecombinationJoules
      + heatOfReactionIonisationJoules;

    translationalEnergy +=
        this->EEle(p) + this->EEle(q) + this->ERot(p)
      + this->ERot(q) + this->EVib(p) + this->EVib(q);

    // centre of mass velocity
    vector Ucm = this->Ucm(p, q);

    label ELevelNewP = 0;
    label ELevelNewQ = 0;
    label vibLevelNewP = 0;
    scalar ERotNewP = 0.0;
    scalar ERotNewQ = 0.0;

    scalar jMaxNewP =
        cloud_.constProps(assIonProductIds[0]).nElectronicLevels();
    scalar jMaxNewQ =
        cloud_.constProps(assIonProductIds[1]).nElectronicLevels();

    const scalarList& EElistNewP =
        cloud_.constProps(assIonProductIds[0]).electronicEnergyList();
    const scalarList& EElistNewQ =
        cloud_.constProps(assIonProductIds[1]).electronicEnergyList();

    const labelList& gListNewP =
        cloud_.constProps(assIonProductIds[0]).degeneracyList();
    const labelList& gListNewQ =
        cloud_.constProps(assIonProductIds[1]).degeneracyList();

    scalar omegaNewPQ =
        0.5
       *(
            cloud_.constProps(assIonProductIds[0]).omega()
          + cloud_.constProps(assIonProductIds[1]).omega()
        );

    ELevelNewP =
        cloud_.postCollisionElectronicEnergyLevel
        (
            translationalEnergy,
            jMaxNewP,
            omegaNewPQ,
            EElistNewP,
            gListNewP
        );

    translationalEnergy -= EElistNewP[ELevelNewP];

    ELevelNewQ =
        cloud_.postCollisionElectronicEnergyLevel
        (
            translationalEnergy,
            jMaxNewQ,
            omegaNewPQ,
            EElistNewQ,
            gListNewQ
        );

    translationalEnergy -= EElistNewQ[ELevelNewQ];

    label rotationalDofNewP =
        cloud_.constProps(assIonProductIds[0]).rotationalDoF();

    if (rotationalDofNewP > VSMALL)
    {
        scalar thetaVNewP =
            cloud_.constProps(assIonProductIds[0]).thetaV()[0];
        scalar thetaDNewP =
            cloud_.constProps(assIonProductIds[0]).thetaD()[0];

        scalar ZrefNewP = cloud_.constProps(assIonProductIds[0]).Zref()[0];
        scalar refTempZvNewP =
            cloud_.constProps(assIonProductIds[0]).TrefZv()[0];

        scalar ChiB = 2.5 - omegaNewPQ;

        if (rotationalDofNewP > VSMALL)
        {
            label iMax =
                translationalEnergy/(physicoChemical::k.value()*thetaVNewP);

            vibLevelNewP =
                cloud_.postCollisionVibrationalEnergyLevel
                (
                    true,
                    0,
                    iMax,
                    thetaVNewP,
                    thetaDNewP,
                    refTempZvNewP,
                    omegaNewPQ,
                    ZrefNewP,
                    translationalEnergy
                );

            translationalEnergy -=
                vibLevelNewP*thetaVNewP*physicoChemical::k.value();

            ERotNewP =
                translationalEnergy
               *cloud_.postCollisionRotationalEnergy(rotationalDofNewP, ChiB);

            translationalEnergy -= ERotNewP;
        }
    }

    scalar mP2 = cloud_.constProps(assIonProductIds[0]).mass();
    scalar mQ2 = cloud_.constProps(assIonProductIds[1]).mass();

    scalar mR = mP2*mQ2/(mP2 + mQ2);

    scalar relVel = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part for collision
    scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    vector postCollisionRelU =
        relVel
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    vector UP = Ucm + (postCollisionRelU*mQ2/(mP2 + mQ2));
    vector UQ = Ucm - (postCollisionRelU*mP2/(mP2 + mQ2));

    vector positionP(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        positionP,
        cell,
        tetFace,
        tetPt
    );

    scalar CWFp = p.CWF();
    scalar RWFp = p.RWF();
    labelList vibLevelP(0, 0);
    if (rotationalDofNewP > VSMALL)
    {
        vibLevelP.resize(1);
        vibLevelP[0] = vibLevelNewP;
    }

    p.position() = positionP;
    p.U() = UP;
    p.CWF() = CWFp;
    p.RWF() = RWFp;
    p.ERot() = ERotNewP;
    p.ELevel() = ELevelNewP;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = assIonProductIds[0];
    p.newParcel() = 0;
    p.vibLevel() = vibLevelP;

    vector positionQ(q.position());

    scalar CWFq = q.CWF();
    scalar RWFq = q.RWF();
    labelList vibLevelQ(0, 0);

    cell = -1;
    tetFace = -1;
    tetPt = -1;

    mesh_.findCellFacePt
    (
        positionQ,
        cell,
        tetFace,
        tetPt
    );

    q.position() = positionQ;
    q.U() = UQ;
    q.CWF() = CWFq;
    q.RWF() = RWFq;
    q.ERot() = ERotNewQ;
    q.ELevel() = ELevelNewQ;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = assIonProductIds[1];
    q.newParcel() = 0;
    q.vibLevel() = vibLevelQ;
}

void Foam::uniGasReaction::exchangePQ
(
    const scalar& heatOfReactionExchJoules,
    const labelPair& exchangeProductIds,
    uniGasParcel& p,
    uniGasParcel& q
)
{
    const label typeIdMol = exchangeProductIds[0];
    const label typeIdAtom = exchangeProductIds[1];

    // change species properties

    const scalar mPExch = cloud_.constProps(typeIdAtom).mass();
    const scalar mQExch = cloud_.constProps(typeIdMol).mass();

    scalar mRExch = mPExch*mQExch/(mPExch + mQExch);

    scalar translationalEnergy = this->translationalEnergy(p, q);

    translationalEnergy +=
        ERot(p) + EVib(p) + EEle(p) + EEle(q)
        + heatOfReactionExchJoules;

    scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

    // Variable Hard Sphere collision part for molecule collision

    scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    vector postCollisionRelU =
        relVelExchMol
        *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

    vector UP =
        Ucm(p, q) + postCollisionRelU*mQExch/(mPExch + mQExch);
    vector UQ =
        Ucm(p, q) - postCollisionRelU*mPExch/(mPExch + mQExch);

    vector positionP(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        positionP,
        cell,
        tetFace,
        tetPt
    );

    scalar CWFp = p.CWF();
    scalar RWFp = p.RWF();
    labelList vibLevelAtom(0, 0);

    p.position() = positionP;
    p.U() = UP;
    p.CWF() = CWFp;
    p.RWF() = RWFp;
    p.ERot() = 0.0;
    p.ELevel() = 0;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = typeIdAtom;
    p.newParcel() = 0;
    p.vibLevel() = vibLevelAtom;

    vector positionQ(q.position());

    cell = -1;
    tetFace = -1;
    tetPt = -1;

    mesh_.findCellFacePt
    (
        positionQ,
        cell,
        tetFace,
        tetPt
    );

    scalar CWFq = q.CWF();
    scalar RWFq = q.RWF();
    labelList vibLevelMol(1, 0);

    q.position() = positionQ;
    q.U() = UQ;
    q.CWF() = CWFq;
    q.RWF() = RWFq;
    q.ERot() = 0.0;
    q.ELevel() = 0;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = typeIdMol;
    q.newParcel() = 0;
    q.vibLevel() = vibLevelMol;
}

void Foam::uniGasReaction::chargeExchangePQ
(
    const scalar& heatOfReactionExchJoules,
    const labelPair& chargeExchangeProductIds,
    uniGasParcel& p,
    uniGasParcel& q
)
{
    scalar translationalEnergy = this->translationalEnergy(p, q);

    translationalEnergy += heatOfReactionExchJoules;

    translationalEnergy +=
        ERot(p) + EVib(p) + EEle(p)
        + ERot(q) + EVib(q) + EEle(q);

    scalar mR =
        cloud_.constProps(chargeExchangeProductIds[0]).mass()
        *cloud_.constProps(chargeExchangeProductIds[1]).mass()
        /(
            cloud_.constProps(chargeExchangeProductIds[0]).mass()
            + cloud_.constProps(chargeExchangeProductIds[1]).mass()
        );

    scalar relVel = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part for collision of
    // molecules
    scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    vector postCollisionRelU =
        relVel
        *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    scalar mP = cloud_.constProps(chargeExchangeProductIds[0]).mass();
    scalar mQ = cloud_.constProps(chargeExchangeProductIds[1]).mass();

    vector UP = Ucm(p, q) + (postCollisionRelU*mQ/(mP + mQ));
    vector UQ = Ucm(p, q) - (postCollisionRelU*mP/(mP + mQ));

    vector positionP(p.position());

    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        positionP,
        cell,
        tetFace,
        tetPt
    );

    scalar CWFp = p.CWF();
    scalar RWFp = p.RWF();
    scalar ERotP = 0.0;
    labelList vibLevelP(0, 0);

    if
    (
        cloud_.constProps(chargeExchangeProductIds[0])
        .vibrationalDoF()
        > VSMALL
    )
    {
        vibLevelP.setSize(1);
        vibLevelP[0] = 0;
    }
    label ELevelP = 0;

    p.position() = positionP;
    p.U() = UP;
    p.CWF() = CWFp;
    p.RWF() = RWFp;
    p.ERot() = ERotP;
    p.ELevel() = ELevelP;
    p.cell() = cell;
    p.tetFace() = tetFace;
    p.tetPt() = tetPt;
    p.typeId() = chargeExchangeProductIds[0];
    p.newParcel() = 0;
    p.vibLevel() = vibLevelP;

    vector positionQ(q.position());

    cell = -1;
    tetFace = -1;
    tetPt = -1;

    mesh_.findCellFacePt
    (
        positionQ,
        cell,
        tetFace,
        tetPt
    );

    scalar CWFq = q.CWF();
    scalar RWFq = q.RWF();
    scalar ERotQ = 0.0;
    labelList vibLevelQ(0, 0);

    if
    (
        cloud_.constProps(chargeExchangeProductIds[1])
        .vibrationalDoF()
        > VSMALL
    )
    {
        vibLevelQ.setSize(1);
        vibLevelQ[0] = 0;
    }
    label ELevelQ = 0;

    q.position() = positionQ;
    q.U() = UQ;
    q.CWF() = CWFq;
    q.RWF() = RWFq;
    q.ERot() = ERotQ;
    q.ELevel() = ELevelQ;
    q.cell() = cell;
    q.tetFace() = tetFace;
    q.tetPt() = tetPt;
    q.typeId() = chargeExchangeProductIds[1];
    q.newParcel() = 0;
    q.vibLevel() = vibLevelQ;
}


bool Foam::uniGasReaction::relax() const
{
    return relax_;
}


const Foam::uniGasCloud& Foam::uniGasReaction::cloud() const
{
    return cloud_;
}


bool Foam::uniGasReaction::outputResults(const label counterIndex)
{
    return writeRatesToTerminal_;
}


// ************************************************************************* //
