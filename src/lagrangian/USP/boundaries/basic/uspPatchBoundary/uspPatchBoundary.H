/*---------------------------------------------------------------------------* \
  =========                 |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::uspPatchBoundary

Description
    Used for applying instantaneous actions on uspParcels colliding with a
    boundary patch.

SourceFiles
    uspPatchBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef uspPatchBoundary_H
#define uspPatchBoundary_H

#include "autoPtr.H"
#include "uspBoundaryBase.H"
#include "uspParcel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class uspCloud;
class Time;

/*---------------------------------------------------------------------------*\
                      Class uspPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class uspPatchBoundary
:
    public uspBoundaryBase
{
protected:

    // Protected data

        //- Faces on the boundary
        labelList faces_;
        scalar patchSurfaceArea_;
        scalar totalPatchSurfaceArea_;

        //- Cells corresponding to the boundary faces
        labelList cells_;

        bool measurePropertiesAtWall_;

        //- Measured fields at the patch
        scalar preIE_;
        vector preIMom_;


public:

    //- Runtime type information
    TypeName("uspPatchBoundary");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uspPatchBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    uspPatchBoundary
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );


    // Selectors

        static autoPtr<uspPatchBoundary> New
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uspPatchBoundary() = default;


    // Member Functions

        virtual void measurePropertiesBeforeControl(const uspParcel& p);
        virtual void measurePropertiesAfterControl
        (
            const uspParcel& p,
            const scalar heatOfReaction
        );

        //- Perform a diffuse wall reflection
        virtual void diffuseReflection
        (
            uspParcel& p,
            const scalar wallTemp,
            const vector& wallVelocity
        );

        //- perform a specular reflection
        virtual void specularReflection(uspParcel& p);

        //- Calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- Control the uspParcel
        virtual void controlParticle
        (
            uspParcel& mol,
            uspParcel::trackingData& td
        ) = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        virtual const labelList& controlPatch() const;

        //- Return the adjacent cells
        virtual const labelList& controlZone() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
