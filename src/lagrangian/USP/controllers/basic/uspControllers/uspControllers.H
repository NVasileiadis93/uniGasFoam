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
    Foam::uspControllers

Description
    Stores all the information for the uspControllers models defined within
    the uspControllersDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef uspControllers_H
#define uspControllers_H

#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "uspStateController.H"
#include "uspFluxController.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class uspControllers Declaration
\*---------------------------------------------------------------------------*/

class uspControllers
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        IOdictionary uspControllersDict_;

        //- Number of state controllers
        label nStateControllers_;

        //- Number of flux controllers
        label nFluxControllers_;

        //- State controllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List<autoPtr<uspStateController>> stateControllers_;

        //- Flux controllers
        PtrList<entry> fluxControllersList_;
        List<word> fCNames_;
        List<label> fCIds_;
        List<fileName> fCFixedPathNames_;
        List<autoPtr<uspFluxController>> fluxControllers_;


public:

    // Constructors

        //- Null Constructor
        uspControllers(const Time& t, const polyMesh& mesh);

        //- Constructor for uspFOAM
        uspControllers(const Time& t, const polyMesh& mesh, uspCloud& cloud);

    //- Destructor
    ~uspControllers() = default;


    // Member Functions

        //- Initial configuration
        void initialConfig();

        //- Update time information
        void updateTimeInfo();

        //- Different control stages
        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        //- Calculate controller properties
        void calculateProps();

        //- Output controller properties
        void outputResults();


        // Access

            label nStateControllers() const;

            inline List<autoPtr<uspStateController>>& stateControllers();

            inline const List<autoPtr<uspStateController>>&
                stateControllers() const;

            inline List<autoPtr<uspFluxController>>& fluxControllers();

            inline const List<autoPtr<uspFluxController>>&
                fluxControllers() const;

            inline const List<word>& stateControllersNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
