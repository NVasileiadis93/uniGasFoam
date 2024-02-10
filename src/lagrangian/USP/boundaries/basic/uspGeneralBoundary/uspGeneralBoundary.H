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
    Foam::uspGeneralBoundary

Description

SourceFiles
    uspGeneralBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef uspGeneralBoundary_H
#define uspGeneralBoundary_H

#include "autoPtr.H"
#include "uspBoundaryBase.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class uspCloud;

/*---------------------------------------------------------------------------*\
                     Class uspGeneralBoundary Declaration
\*---------------------------------------------------------------------------*/

class uspGeneralBoundary
:
    public uspBoundaryBase
{
protected:

    // Protected data

        //- Faces on the boundary
        labelList faces_;

        //- Type IDs (species)
        labelList typeIds_;

        //- Patch surface area
        scalar patchSurfaceArea_;

        //- Cells corresponding to the boundary faces
        labelList cells_;

        // Set all the properties below from model if required
        List<scalarField> accumulatedParcelsToInsert_;


public:

    //- Runtime type information
    TypeName("uspGeneralBoundary");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uspGeneralBoundary,
        dictionary,
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        ),
        (mesh, cloud, dict)
    );

    //- Construct from components
    uspGeneralBoundary
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );


    //- Selectors
    static autoPtr<uspGeneralBoundary> New
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~uspGeneralBoundary() = default;


    // Member Functions

        //- Create an initial configuration
        virtual void initialConfiguration() = 0;

        //- Control
        virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- Calculate any required properties
        virtual void calculateProperties() = 0;

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalarList& numDen,
            const scalar& transT,
            const vector& velocity
        );

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalarList& numDen,
            const scalar& transT,
            const vector& velocity,
            const vector& heatFlux,
            const tensor& stress
        );

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const List<scalarField>& numDen,
            const scalarField& transT,
            const vectorField& velocity
        );

        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const List<scalarField>& numDen,
            const scalarField& transT,
            const vectorField& velocity,
            const vectorField& heatFlux,
            const tensorField& stress
        );
               
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalar& numDen,
            const scalarField& molFractions,
            const scalar& transT,
            const vectorField& velocity
        );
        
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const scalarField& numDen,
            const scalarField& molFractions,
            const scalarField& transT,
            const vectorField& velocity
        );
      
        //- Calculate how many parcels to insert
        virtual void computeParcelsToInsert
        (
            const List<scalarField>& numDen,
            const scalar& transT,
            const vectorField& velocity
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& transT,
            const scalar& rotT,
            const scalar& vibT,
            const scalar& elecT,
            const vector& velocity
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalarList& numDen,
            const scalar& transT,
            const scalar& rotT,
            const scalar& vibT,
            const scalar& elecT,
            const vector& velocity,
            const vector& heatFlux,
            const tensor& stress
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalarField& transT,
            const scalarField& rotT,
            const scalarField& vibT,
            const scalarField& elecT,
            const vectorField& velocity
        );
 
        //- Insert parcels
        virtual void insertParcels
        (
            const List<scalarField>& numDen,
            const scalarField& transT,
            const scalarField& rotT,
            const scalarField& vibT,
            const scalarField& elecT,
            const vectorField& velocity,
            const vectorField& heatFlux,
            const tensorField& stress
        );

        //- Insert parcels
        virtual void insertParcels
        (
            const scalar& transT,
            const vectorField& velocity
        );
        
        //- Insert parcels
        virtual void insertParcels
        (
            const scalarField& transT,
            const vectorField& velocity
        );

        // Access

            virtual const labelList& controlPatch() const;

            //- Return the adjacent cells
            virtual const labelList& controlZone() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
