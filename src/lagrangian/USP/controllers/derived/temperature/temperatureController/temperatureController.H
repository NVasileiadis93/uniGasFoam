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
    temperatureController

Description
    Controls translational temperature on a cell-by-cell basis in a user-defined
    zone on a mesh. Provides the best resolution, but the worst statistics (low
    number of particles per cell). No need for parallel processing as everything
    is on a per-cell basis.

SourceFiles
    temperatureController.C


\*---------------------------------------------------------------------------*/

#ifndef temperatureController_H
#define temperatureController_H

#include "uspStateController.H"
#include "vector.H"
#include "vectorField.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class temperatureController Declaration
\*---------------------------------------------------------------------------*/

class temperatureController
:
   public uspStateController
{
    // Private Data

        //- subDictionary containing the properties
        dictionary propsDict_;

        dictionary timeDictVel_;
        timeDataMeas timeVel_;

        //- time constant
        scalar tauT_;

        // component control

        bool componentControl_;
        bool X_;
        bool Y_;
        bool Z_;

        labelList typeIds_;

        //- measurements
        scalarField massV_;
        vectorField momV_;
        vectorField UMean_;
        scalarField mcc_;
        scalarField m_;
        scalarField nParcels_;

        scalarField measuredTranslationalTemperature_;

        //- velocity rescaling term
        scalarField chi_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("temperatureController");

    // Constructors

        //- Construct from components
        temperatureController
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

        //- Destructor
        virtual ~temperatureController() = default;


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
        virtual void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
