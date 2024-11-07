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
    Foam::uniGasFaceTracker

Description
    Class for tracking particles microscopic properties while crossing mesh faces.
    This is used to calculate particle, mass and momentum flow rate.

SourceFiles
    uniGasFaceTrackerI.H
    uniGasFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasFaceTracker_H
#define uniGasFaceTracker_H

#include "vector.H"
#include "fvMesh.H"
#include "uniGasParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class uniGasFaceTracker Declaration
\*---------------------------------------------------------------------------*/

class uniGasFaceTracker
{
    // Private data

        const fvMesh& mesh_;

        uniGasCloud& cloud_;

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;
        List<vectorField> momentumIdFlux_;
        List<scalarField> energyIdFlux_;

public:

    friend class uniGasParcel;

    // Constructors

        //- Construct from mesh and cloud
        uniGasFaceTracker
        (
            const polyMesh&,
            uniGasCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (uniGasFoam)
        uniGasFaceTracker
        (
            const polyMesh&,
            uniGasCloud& cloud,
            const bool dummy
        );

    //- Destructor
    ~uniGasFaceTracker() = default;


    // Member Functions

        void updateFields(uniGasParcel& p);

        void clean();


        // Access

            inline const List<scalarField>& parcelIdFlux() const;

            inline const List<scalarField>& massIdFlux() const;

            inline const List<vectorField>& momentumIdFlux() const;

            inline const List<scalarField>& energyIdFlux() const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniGasFaceTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
