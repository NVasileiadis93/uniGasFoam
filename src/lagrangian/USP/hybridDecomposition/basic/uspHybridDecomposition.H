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
    Foam::uspHybridDecomposition

Description
    Templated hybrid decomposition class

SourceFiles
    uspHybridDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef uspHybridDecomposition_H
#define uspHybridDecomposition_H

#include "fvMesh.H"
#include "uspParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class uspCloud;
class uspParcel;

/*---------------------------------------------------------------------------*\
                    Class uspHybridDecomposition Declaration
\*---------------------------------------------------------------------------*/

class uspHybridDecomposition
{
    // Private data

protected:

        //- reference to time
        const Time& time_;

        //- reference to the mesh
        const fvMesh& mesh_;

        // reference to the owner cloud class
        uspCloud& cloud_;

        //- The hybrid decomposition IO dictionary
        IOdictionary hybridDecompositionDict_;

        //- The time properties dictionary
        dictionary timeDict_;

        // Time interval between successive decompositions
        label decompositionInterval_;
        
        // Breakdown parameter threshold
        bool resetAtDecomposition_;
        
        // Time averaging coefficient
        scalar resetAtDecompositionUntilTime_;

        // Number of refinement passes
        label refinementPasses_;

        // Levels of neightbors
        label neighborLevels_;

        // Fraction of maximum allowed neighbors of opposite collision model
        scalar maxNeighborFraction_;

        // Coefficient for correcting boundary measurements
        scalar boundCoeff_;

        // List of neighborhood cells
        DynamicList<label> neighborCells_;

public:


    //- Default name of the hybrid decomposition dictionary
    static const word dictName;

    //- Runtime type information
    TypeName("uspHybridDecomposition");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uspHybridDecomposition,
        dictionary,
        (
            const Time& t,
            const polyMesh& mesh,
            uspCloud& owner
        ),
        (t, mesh, owner)
    );

    // Constructors

        //- Construct from components
        uspHybridDecomposition
        (
            const Time& t,
            const polyMesh& mesh,
            uspCloud& owner
        );

    //- Destructor
    virtual ~uspHybridDecomposition() = default;


    //- Selector
    static autoPtr<uspHybridDecomposition> New
    (
        const Time& t,
        const polyMesh& mesh,
        uspCloud& owner
    );

    // Access

    // Member Functions

        void fetchCellNeighborhood
        (
            const label cell,
            const label nLevels,
            DynamicList<label>& neighborCells
        );

        void updateProperties();

        virtual void decompose() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
