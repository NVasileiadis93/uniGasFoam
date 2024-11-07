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

#include "uniGasField.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniGasField, 0);

defineRunTimeSelectionTable(uniGasField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniGasField::uniGasField
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    timeDict_(dict.subDict("timeProperties")),
    sampleInterval_(timeDict_.getOrDefault<label>("sampleInterval",1)),
    resetFieldsAtOutput_(timeDict_.getOrDefault<bool>("resetAtOutput",true)),
    resetFieldsAtOutputUntilTime_(timeDict_.getOrDefault<scalar>("resetAtOutputUntilTime",VGREAT)),
    casePath_(t.path()/"fieldMeasurements"),
    timePath_()
{
    // Needed?
    mkDir(casePath_);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<uniGasField> uniGasField::New
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("fieldModel"));

    Info<< "Selecting field: " << modelType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uniGasField",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uniGasField>(cstrIter()(t, mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uniGasField::updateProperties(const dictionary& dict)
{

    timeDict_ = dict.subDict("timeProperties");

    timeDict_.readIfPresent("sampleInterval", sampleInterval_);

    timeDict_.readIfPresent("resetAtOutput", resetFieldsAtOutput_);

    timeDict_.readIfPresent("resetAtOutputUntilTime", resetFieldsAtOutputUntilTime_);

}

const label& uniGasField::sampleInterval() const
{
    return sampleInterval_;
}


const bool& uniGasField::resetFieldsAtOutput() const
{
    return resetFieldsAtOutput_;
}


const scalar& uniGasField::resetFieldsAtOutputUntilTime() const
{
    return resetFieldsAtOutputUntilTime_;
}


const fileName& uniGasField::casePath() const
{
    return casePath_;
}


fileName& uniGasField::casePath()
{
    return casePath_;
}


const fileName& uniGasField::timePath() const
{
    return timePath_;
}


fileName& uniGasField::timePath()
{
    return timePath_;
}

} // End namespace Foam

// ************************************************************************* //
