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

#include "uniGasAllConfigurations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::uniGasAllConfigurations::uniGasAllConfigurations
(
    const IOdictionary& uniGasInitialisationDict,
    uniGasCloud& cloud
)
:
    uniGasInitialisationDict_(uniGasInitialisationDict),
    configurationList_(uniGasInitialisationDict_.lookup("configurations")),
    ids_(configurationList_.size()),
    configurations_(configurationList_.size())
{
    Info << nl << "Creating uniGas configurations: " << nl << endl;

    if (configurations_.size() > 0)
    {
        forAll(configurations_, c)
        {
            const entry& configurationI = configurationList_[c];
            const dictionary& configurationIDict = configurationI.dict();

            configurations_[c] =
                uniGasConfiguration::New(cloud, configurationIDict);

            ids_[c] = c;
        }
    }
}


void Foam::uniGasAllConfigurations::setInitialConfig()
{
    forAll(configurations_, c)
    {
        configurations_[c]->setInitialConfiguration();
    }
}


// ************************************************************************* //
