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

Class
    Foam::uspAllConfigurations

Description

\*---------------------------------------------------------------------------*/

#ifndef uspAllConfigurations_H
#define uspAllConfigurations_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "uspConfiguration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class uspAllConfigurations Declaration
\*---------------------------------------------------------------------------*/

class uspAllConfigurations
{
    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        const IOdictionary& uspInitialiseDict_;

        PtrList<entry> configurationList_;
        List<label> ids_;
        List<autoPtr<uspConfiguration>> configurations_;


public:

    // Constructors

        //- Null Constructor
        uspAllConfigurations
        (
            const IOdictionary& uspInitialiseDict,
            uspCloud& cloud
        );


    // Member Functions

        //- initial configuration
        void setInitialConfig();

        // Access

            inline List<autoPtr<uspConfiguration>>& configurations();
            inline const List<autoPtr<uspConfiguration>>&
                configurations() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspAllConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
