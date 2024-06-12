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
    Foam::uniGasPatchBoundary

Description
    Used for applying instantaneous actions on uniGasParcels colliding with a
    boundary patch.

SourceFiles
    uniGasPatchBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasPatchBoundary_H
#define uniGasPatchBoundary_H

#include "autoPtr.H"
#include "uniGasBoundaryBase.H"
#include "uniGasParcel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class uniGasCloud;
class Time;

/*---------------------------------------------------------------------------*\
                      Class uniGasPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class uniGasPatchBoundary
:
    public uniGasBoundaryBase
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
    TypeName("uniGasPatchBoundary");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uniGasPatchBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    uniGasPatchBoundary
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );


    // Selectors

        static autoPtr<uniGasPatchBoundary> New
        (
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uniGasPatchBoundary() = default;


    // Member Functions

        virtual void measurePropertiesBeforeControl(const uniGasParcel& p);
        virtual void measurePropertiesAfterControl
        (
            const uniGasParcel& p,
            const scalar heatOfReaction
        );

        //- Perform a diffuse wall reflection
        virtual void diffuseReflection
        (
            uniGasParcel& p,
            const scalar wallTemp,
            const vector& wallVelocity
        );

        //- perform a specular reflection
        virtual void specularReflection(uniGasParcel& p);

        //- Calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- Control the uniGasParcel
        virtual void controlParticle
        (
            uniGasParcel& mol,
            uniGasParcel::trackingData& td
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
