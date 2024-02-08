/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    uspMassFluxSurface

Description
    measure the average mass flux at a user-defined surface

Note
    To understand this measurement model, one needs to look at how the cloudTracker works,
    as well as how processor patches are defined on the mesh.

    - cloud tracker fields are generated every time-step (used by measurement models)
      then fields are cleaned.

    - in parallel processing runs a facezone on a processor patch does not get
      divided uniquely. A pair of coupled processor patches will both identify that
      a face belonging to the face occurs. care must therefore be taken to get the flux
      right. For example note how the surface area of the total face zone is calculated.
      Faces on processor patches are given a weighting of 0.5 so that when communication
      occurs the correct zonal surface area is computed. The same is said for atomisticMolecule
      crossings. In the cloudTracker class only one face on the processor patch gets
      gets assigned the flux. Thereby synchronisation between processors gives the
      correct total flux.

SourceFiles
    uspMassFluxSurface.C

\*---------------------------------------------------------------------------*/

#ifndef uspMassFluxSurface_H
#define uspMassFluxSurface_H

#include "uspField.H"
#include "vector.H"
#include "uspCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspMassFluxSurface Declaration
\*---------------------------------------------------------------------------*/

class uspMassFluxSurface
:
	public uspField
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;

        word fieldName_;

        label regionId_;
        label sampleInterval_;
        label sampleCounter_;
        word faceZoneName_;

        scalar zoneSurfaceArea_;

        //- Parcel IDs
        labelList typeIds_;

        vector fluxDirection_;

        scalar molsZone_;
        scalar massZone_;
        scalar momentumZone_;
        scalar averagingCounter_;
        label timeIndex_;

        scalarField molFluxZone_;
        scalarField massFluxZone_;
        scalarField massFlowZone_;
        scalarField momentumFlowZone_;

        bool averagingAcrossManyRuns_;


    // Private Member Functions

        void read();

        void write();


 public:

    //- Runtime type information
    TypeName("uspMassFluxSurface");

    // Constructors

        //- Construct from components
        uspMassFluxSurface
        (
            const Time& t,
            const polyMesh& mesh,
            uspCloud& molCloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uspMassFluxSurface() = default;


    // Member Functions

        //- Create initial field
        virtual void createField();

        //- Calculate field
        virtual void calculateField();


        //- Write field
        virtual void writeField();

        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
