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
    Foam::densityController

Description
    This is a density controller, which uses a slightly modified version of the
    particle insertion protocol called USHER, proposed by
    Delgado-Buscalioni and Coveney.

SourceFiles
    densityController.C

\*---------------------------------------------------------------------------*/

#ifndef densityController_H
#define densityController_H

#include "uspStateController.H"
#include "vector.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class densityController Declaration
\*---------------------------------------------------------------------------*/

class densityController
:
    public uspStateController
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;
        label typeId_;

        //- Initial target density, velocity and temperature
        scalar density_; // parcel density
        vector velocity_;
        scalar temperature_;
        scalarField densities_;
        vectorField velocities_;
        scalarField temperatures_;
        labelField nMolsToControl_;
        scalarField parcels_;
        scalarField measuredDensity_;
        scalarField residual_;
        scalarField residualSum_;


    // Private Member Functions

        void deleteParcelFromUSP(const label c);
        void insertParcelWithinUSP(const label c);

        void setProperties();


public:

    //- Runtime type information
    TypeName("densityController");

    //- Construct from components
    densityController
    (
        const polyMesh& mesh,
        uspCloud& uspCloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~densityController() = default;


    // Member Functions

        virtual void initialConfiguration();

        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
