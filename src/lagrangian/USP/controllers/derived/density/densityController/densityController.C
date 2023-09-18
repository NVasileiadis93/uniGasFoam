/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "densityController.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(densityController, 0);
addToRunTimeSelectionTable(uspStateController, densityController, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityController::densityController
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspStateController(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeId_(-1),
    density_(propsDict_.get<scalar>("massDensity")),
    velocity_(propsDict_.get<vector>("velocity")),
    temperature_(propsDict_.get<scalar>("temperature")),
    densities_(controlZone().size(), Zero),
    velocities_(controlZone().size(), Zero),
    temperatures_(controlZone().size(), Zero),
    nMolsToControl_(controlZone().size(), Zero),
    parcels_(controlZone().size(), Zero),
    measuredDensity_(controlZone().size(), Zero),
    residual_(controlZone().size(), Zero),
    residualSum_(controlZone().size(), Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::densityController::initialConfiguration()
{
    setProperties();
}


void Foam::densityController::calculateProperties()
{
    if (timeData_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        label c = 0;
        for (const label celli : controlZone())
        {
            for (uspParcel* p : cellOccupancy[celli])
            {
                if (p->typeId() == typeId_)
                {
                    scalar CWF =
                        cloud_.cellWF(celli);
                    scalar RWF =
                        cloud_.axiRWF(cloud_.mesh().cellCentres()[celli]);

                    parcels_[c] += CWF*RWF;
                }
            }
            ++c;
        }
    }

    if (timeData_.averagingTime() && (controlZone().size() > 0))
    {
        Info<< "Averaging" << endl;
        const scalar nAvTimeSteps = timeData_.nAvTimeSteps().value();
        Info<< "nAvTimeSteps = " << nAvTimeSteps << endl;

        forAll(measuredDensity_, c)
        {
            measuredDensity_[c] = parcels_[c]/nAvTimeSteps;
            Info<< "measuredDensity_[c] = " << measuredDensity_[c] << endl;
            Info<< "densities_[c] = " << densities_[c] << endl;
            scalar deltaNParcels = densities_[c] - measuredDensity_[c];
            Info<< "deltaNParcels = " << deltaNParcels << endl;
            scalar integerNParcels = 0.0;

            if (deltaNParcels > 0.0)
            {
                integerNParcels = label(deltaNParcels + 0.5);
            }
            else if (deltaNParcels < 0.0)
            {
                integerNParcels = label(deltaNParcels - 0.5);
            }

            nMolsToControl_[c] = integerNParcels;
            residual_[c] = deltaNParcels - nMolsToControl_[c];
            parcels_[c] = 0.0;

            Info<< "Number of parcels to control = "
                << integerNParcels << nl
                << "Number of residual parcels = " << residual_[c] << endl;
        }
    }
}


void Foam::densityController::controlParcelsBeforeMove()
{
    if (control_ && timeData_.controlTime())
    {
        Info<< "Controlling density" << endl;

        labelField nMolsToControl(nMolsToControl_.size(), 0);
        scalarField residualSum(residualSum_.size(), 0);

        forAll(nMolsToControl_, c)
        {
            if (nMolsToControl_[c] > 0) // insert parcels
            {
                insertParcelWithinUSP(c);
            }
            else if (nMolsToControl_[c] < 0) // delete parcels
            {
                deleteParcelFromUSP(c);
            }
        }

        forAll(residualSum_, c)
        {
            residualSum_[c] += residual_[c];

            if (residualSum_[c] > 1) // insert residual parcel
            {
                insertParcelWithinUSP(c);
                residualSum_[c] -= 1;
            }
            else if (residualSum_[c] < -1) // delete residual parcel
            {
                deleteParcelFromUSP(c);
                residualSum_[c] += 1;
            }
        }

        nMolsToControl_ = 0;
    }
}


void Foam::densityController::deleteParcelFromUSP(const label c)
{
    const label cellI = controlZone()[c];
    const auto& cellOccupancy = cloud_.cellOccupancy();
    const auto& molsInCell = cellOccupancy[cellI];

    if (molsInCell.size() > 0)
    {
        label removei = rndGen_.position<label>(0, molsInCell.size() - 1);
        uspParcel* delParcel = molsInCell[removei];

        // delete molecule from cellOccupancy (before deleting it from cloud)
        cloud_.removeParcelFromCellOccupancy(removei, cellI);
        cloud_.deleteParticle(*delParcel);
    }
}


void Foam::densityController::insertParcelWithinUSP(const label c)
{
    const label cellI = controlZone()[c];

    const vector& cC = mesh_.cellCentres()[cellI];

    // find the maximum distance between cell centre and cell vertices
    const labelList& cellPoints = mesh_.cellPoints()[cellI];
    scalar maxDistance = 0.0;

    forAll(cellPoints, cP)
    {
        const vector& vertexI = mesh_.points()[cellPoints[cP]];
        maxDistance = max(maxDistance, mag(vertexI - cC));
    }

    // find a random point within the cell
    bool isPointInCell = false;

    vector p(Zero);

    while (!isPointInCell)
    {
        // select and normalise random direction
        vector randDirection(rndGen_.GaussNormal<vector>());
        randDirection.normalise();

        p = randDirection*rndGen_.sample01<scalar>()*maxDistance + cC;

        if (mesh_.pointInCell(p, cellI))
        {
            isPointInCell = true;
        }
    }

    const uspParcel::constantProperties& cP = cloud_.constProps(typeId_);

    vector U = cloud_.equipartitionLinearVelocity(temperature_, cP.mass());

    scalar ERot = cloud_.equipartitionRotationalEnergy
    (
        temperature_,
        cP.rotationalDoF()
    );

    labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
    (
        temperature_,
        cP.vibrationalDoF(),
        typeId_
    );

    label ELevel = cloud_.equipartitionElectronicLevel
    (
        temperature_,
        cP.degeneracyList(),
        cP.electronicEnergyList(),
        typeId_
    );

    // thermal velocity + stream velocity = instantaneous velocity -
    // STREAM VELOCITY HERE IS USER DEFINED in the controllersDict.
    U += velocity_;


    scalar CWF = cloud_.cellWF(cellI);
    scalar RWF = cloud_.axiRWF(mesh_.cellCentres()[cellI]);

    cloud_.addNewParcel
    (
        p,
        U,
        CWF,
        RWF,
        ERot,
        ELevel,
        cellI,
        typeId_,
        0,
        vibLevel
    );
}


void Foam::densityController::controlParcelsBeforeCollisions()
{}


void Foam::densityController::controlParcelsAfterCollisions()
{}


void Foam::densityController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::densityController::updateProperties(const dictionary& dict)
{
    // the main controller properties should be updated first
    uspStateController::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    if (readStateFromFile_)
    {
        density_ = propsDict_.get<scalar>("massDensity");
        velocity_ = propsDict_.get<vector>("velocity");
        temperature_ = propsDict_.get<scalar>("temperature");

        density_ /= cloud_.constProps(typeId_).mass()*cloud_.nParticle();

        const scalarField& volField = mesh_.cellVolumes();

        forAll(densities_, c)
        {
            const label cellI = controlZone()[c];
            densities_[c] = density_*volField[cellI];
        }

        forAll(velocities_, c)
        {
            velocities_[c] = velocity_;
        }

        forAll(temperatures_, c)
        {
            temperatures_[c] = temperature_;
        }
    }
}


void Foam::densityController::setProperties()
{
    // search for mol id
    const word typeId = propsDict_.get<word>("typeId");
    typeId_ = cloud_.typeIdList().find(typeId);

    if (typeId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find typeId: " << typeId << nl << "in: "
            << mesh_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    Info<< "typeId_ = " << typeId_ << endl;
    density_ /= cloud_.constProps(typeId_).mass()*cloud_.nParticle();
    const scalarField& volField = mesh_.cellVolumes();

    forAll(densities_, c)
    {
        const label cellI = controlZone()[c];
        densities_[c] = density_*volField[cellI];
        Info<< "densities_[c] = " << densities_[c] << endl;
    }

    forAll(velocities_, c)
    {
        velocities_[c] = velocity_;
    }

    forAll(temperatures_, c)
    {
        temperatures_[c] = temperature_;
    }
}


// ************************************************************************* //
