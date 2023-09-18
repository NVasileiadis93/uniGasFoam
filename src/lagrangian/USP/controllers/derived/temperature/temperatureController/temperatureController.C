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

#include "temperatureController.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(temperatureController, 0);

addToRunTimeSelectionTable
(
    uspStateController,
    temperatureController,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureController::temperatureController
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspStateController(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeDictVel_(propsDict_.subDict("timePropertiesForVelocity")),
    timeVel_(mesh_.time(), timeDictVel_),
    tauT_(propsDict_.get<scalar>("tauT")),
    componentControl_(false),
    X_(false),
    Y_(false),
    Z_(false),
    typeIds_(),
    massV_(controlZone().size(), Zero),
    momV_(controlZone().size(), Zero),
    UMean_(controlZone().size(), Zero),
    mcc_(controlZone().size(), Zero),
    m_(controlZone().size(), Zero),
    nParcels_(controlZone().size(), Zero),
    measuredTranslationalTemperature_(controlZone().size(), Zero),
    chi_(controlZone().size(), Zero)
{
    if (propsDict_.readIfPresent("componentControl", componentControl_))
    {
        if (componentControl_)
        {
            X_ = propsDict_.found("X");
            Y_ = propsDict_.found("Y");
            Z_ = propsDict_.found("Z");

            Info << "X_ = " << X_ << ", Y_ = " << Y_ << ", Z = " << Z_ << endl;

            if (!X_ && !Y_ && !Z_)
            {
                FatalErrorInFunction
                    << "At least one component (X, Y, Z) should be chosen "
                    << nl << "in: "
                    << mesh_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }

    setProperties();

    measuredTranslationalTemperature_ = temperature_;

    typeIds_ = cloud_.getTypeIDs(propsDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperatureController::initialConfiguration()
{}


void Foam::temperatureController::calculateProperties()
{
    timeVel_++;

    const labelList& cells = mesh_.cellZones()[regionId_];

    if (timeVel_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        forAll(controlZone(), c)
        {
            const label cellI = cells[c];

            for (uspParcel* p : cellOccupancy[cellI])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    scalar CWF =
                        cloud_.cellWF(cellI);
                    scalar RWF =
                        cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);

                    scalar nParticle = cloud_.nParticle()*CWF*RWF;

                    const scalar mass =
                        cloud_.constProps(p->typeId()).mass()*nParticle;

                    massV_[c] += mass;
                    momV_[c] += p->U()*mass;
                }
            }
        }
    }

    if (timeVel_.averagingTime())
    {
        UMean_ = vector::zero;

        forAll(UMean_, c)
        {
            if (massV_[c] > 0.0)
            {
                UMean_[c] = momV_[c]/massV_[c];
            }
        }

        // reset
        if (timeData_.resetFieldsAtOutput())
        {
            massV_ = 0.0;
            momV_ = vector::zero;
        }
    }

    if (timeData_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        forAll(controlZone(), c)
        {
            const label cellI = cells[c];

            for (uspParcel* p : cellOccupancy[cellI])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    scalar CWF =
                        cloud_.cellWF(cellI);
                    scalar RWF =
                        cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                    scalar nParticle = cloud_.nParticle()*CWF*RWF;

                    scalar mass =
                        cloud_.constProps(p->typeId()).mass()*nParticle;

                    mcc_[c] += mass*magSqr(p->U());
                    m_[c] += mass;
                    nParcels_[c] += nParticle;
                }
            }
        }
    }

    if (timeData_.averagingTime())
    {
        measuredTranslationalTemperature_ = scalar(0);

        const scalar deltaTUSP = mesh_.time().deltaTValue();

        forAll(measuredTranslationalTemperature_, c)
        {
            if (nParcels_[c] > 0)
            {
                measuredTranslationalTemperature_[c] =
                    (1.0/(3.0*physicoChemical::k.value()))
                   *(
                        (mcc_[c]/nParcels_[c])
                      - ((m_[c]/nParcels_[c])*magSqr(UMean_[c]))
                    );

                chi_[c] =
                    sqrt
                    (
                        1.0
                      + (deltaTUSP/tauT_)
                       *(
                            temperature_/measuredTranslationalTemperature_[c]
                          - 1.0
                        )
                    );

                Info<< "target temperature: " << temperature_
                    << " UMean_ : " << UMean_[c]
                    << " measured T: " << measuredTranslationalTemperature_[c]
                    << " chi: " << chi_[c]
                    << endl;
            }
        }

        // Reset
        if (timeData_.resetFieldsAtOutput())
        {
            mcc_ = 0.0;
            m_ = 0.0;
            nParcels_ = 0.0;
        }
    }
}


void Foam::temperatureController::controlParcelsBeforeMove()
{
    if (control_ && timeData_.controlTime())
    {
        Info << "temperatureController: control" << endl;

        const labelList& cells = mesh_.cellZones()[regionId_];

        const auto& cellOccupancy = cloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label celli = cells[c];

            for (uspParcel* p : cellOccupancy[celli])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    if (componentControl_)
                    {
                        if (X_ && chi_[c] > 0)
                        {
                            p->U().x() *= chi_[c];
                        }
                        if (Y_ && chi_[c] > 0)
                        {
                            p->U().y() *= chi_[c];
                        }
                        if (Z_ && chi_[c] > 0)
                        {
                            p->U().z() *= chi_[c];
                        }
                    }
                    else
                    {
                        if (chi_[c] > 0)
                        {
                            p->U() -= UMean_[c];
                            p->U() *= chi_[c];
                            p->U() += UMean_[c];
                        }
                    }
                }
            }
        }
    }
}


void Foam::temperatureController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::temperatureController::controlParcelsBeforeCollisions()
{}


void Foam::temperatureController::controlParcelsAfterCollisions()
{}


void Foam::temperatureController::updateProperties(const dictionary& dict)
{
    // The main controller properties should be updated first
    uspStateController::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");
}


void Foam::temperatureController::setProperties()
{
    temperature_ = propsDict_.get<scalar>("controlTemperature");
}


// ************************************************************************* //
