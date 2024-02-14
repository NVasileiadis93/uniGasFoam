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

#include "uspField.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uspField, 0);

defineRunTimeSelectionTable(uspField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uspField::uspField
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud,
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

autoPtr<uspField> uspField::New
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud,
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
            "uspField",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspField>(cstrIter()(t, mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uspField::updateProperties(const dictionary& dict)
{

    timeDict_ = dict.subDict("timeProperties");

    timeDict_.readIfPresent("sampleInterval", sampleInterval_);

    timeDict_.readIfPresent("resetAtOutput", resetFieldsAtOutput_);

    timeDict_.readIfPresent("resetAtOutputUntilTime", resetFieldsAtOutputUntilTime_);

}

const label& uspField::sampleInterval() const
{
    return sampleInterval_;
}


const bool& uspField::resetFieldsAtOutput() const
{
    return resetFieldsAtOutput_;
}


const scalar& uspField::resetFieldsAtOutputUntilTime() const
{
    return resetFieldsAtOutputUntilTime_;
}


const fileName& uspField::casePath() const
{
    return casePath_;
}


fileName& uspField::casePath()
{
    return casePath_;
}


const fileName& uspField::timePath() const
{
    return timePath_;
}


fileName& uspField::timePath()
{
    return timePath_;
}

} // End namespace Foam

// ************************************************************************* //
