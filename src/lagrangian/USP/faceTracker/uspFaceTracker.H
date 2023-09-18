/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::uspFaceTracker

Description

SourceFiles
    uspFaceTrackerI.H
    uspFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef uspFaceTracker_H
#define uspFaceTracker_H

#include "vector.H"
#include "fvMesh.H"
#include "uspParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class uspFaceTracker Declaration
\*---------------------------------------------------------------------------*/

class uspFaceTracker
{
    // Private data

        const fvMesh& mesh_;

        uspCloud& cloud_;

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;
        List<vectorField> momentumIdFlux_;

public:

    friend class uspParcel;

    // Constructors

        //- Construct from mesh and cloud
        uspFaceTracker
        (
            const polyMesh&,
            uspCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (uspFoam)
        uspFaceTracker
        (
            const polyMesh&,
            uspCloud& cloud,
            const bool dummy
        );

    //- Destructor
    ~uspFaceTracker() = default;


    // Member Functions

        void updateFields(uspParcel& p);

        void clean();


        // Access

            inline const List<scalarField>& parcelIdFlux() const;

            inline const List<scalarField>& massIdFlux() const;

            inline const List<vectorField>& momentumIdFlux() const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspFaceTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
