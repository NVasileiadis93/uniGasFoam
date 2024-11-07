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

Class
    Foam::uniGasAllConfigurations

Description
    Templated class for initialising particles.

Source Files
    uniGasAllConfigurationsI.H
    uniGasAllConfigurations.C
\*---------------------------------------------------------------------------*/

#ifndef uniGasAllConfigurations_H
#define uniGasAllConfigurations_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "uniGasConfiguration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class uniGasAllConfigurations Declaration
\*---------------------------------------------------------------------------*/

class uniGasAllConfigurations
{
    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        const IOdictionary& uniGasInitialisationDict_;

        PtrList<entry> configurationList_;
        List<label> ids_;
        List<autoPtr<uniGasConfiguration>> configurations_;


public:

    // Constructors

        //- Null Constructor
        uniGasAllConfigurations
        (
            const IOdictionary& uniGasInitialisationDict,
            uniGasCloud& cloud
        );


    // Member Functions

        //- initial configuration
        void setInitialConfig();

        // Access

            inline List<autoPtr<uniGasConfiguration>>& configurations();
            inline const List<autoPtr<uniGasConfiguration>>&
                configurations() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniGasAllConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
