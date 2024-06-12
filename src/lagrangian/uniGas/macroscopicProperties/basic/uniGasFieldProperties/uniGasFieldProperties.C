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

#include "uniGasFieldProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasFieldProperties::uniGasFieldProperties
(
    const Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    uniGasFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(),
    fieldNames_(),
    fieldIds_(),
    fields_()
{}


Foam::uniGasFieldProperties::uniGasFieldProperties
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud
)
:
    time_(t),
    uniGasFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(uniGasFieldPropertiesDict_.lookup("uniGasFields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size())
{
    if (fields_.size() > 0)
    {
        Info << "Creating fields: " << nl << endl;

        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();

            fields_[f] = uniGasField::New(time_, mesh, cloud, fieldIDict);

            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::uniGasFieldProperties::createFields()
{
    Info << nl << "Initialising the measurement fields" << nl << endl;

    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}


void Foam::uniGasFieldProperties::calculateFields()
{
    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}


void Foam::uniGasFieldProperties::writeFields()
{
    // Note - not all fields automatically write out to file
    const Time& runTime = time_;

    fileName timePath(runTime.path()/runTime.timeName()/"uniform");

    if (runTime.writeTime())
    {
        if (Pstream::master())
        {
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }
        }
    }

    forAll(fields_, i)
    {
        fields_[i]->timePath() = timePath;
        fields_[i]->writeField();
    }

    if (runTime.writeTime())
    {
        // Checking for modifications in the IOdictionary
        // this allows for run - time tuning of any parameters.

        // NOTES:
        // At the moment, the dictionary is forced to be re-read every
        // write-interval, and properties within the abstract and models are
        // re-set to zero.
        // The "ideal" case is to have the code identify when the dictionary has
        // been modified, before re-reading it in again. Unfortunately the
        // .modified() function is not working properly.

        fieldList_.clear();

        fieldList_ = uniGasFieldPropertiesDict_.lookup("uniGasFields");

        forAll(fields_, i)
        {
            const entry& fieldEntry = fieldList_[i];
            fields_[i]->updateProperties(fieldEntry.dict());
        }
    }
}


// ************************************************************************* //
