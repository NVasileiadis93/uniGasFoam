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

        //- The cloud dictionary
        const dictionary& dict_;

        //- reference to the mesh
        const fvMesh& mesh_;

        // reference to the owner cloud class
        uspCloud& cloud_;

        // Time interval between successive decompositions
        label decompositionInterval_;
        
        // Breakdown parameter threshold
        scalar breakdownMax_;
        
        // Time averaging coefficient
        scalar theta_;
        
        // Number of smoothing passes
        label smoothingPasses_;

        // Number of refinement passes
        label refinementPasses_;

        // Coefficient for correcting boundary measurements
        scalar boundCoeff_;


public:

    //- Runtime type information
    TypeName("uspHybridDecomposition");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uspHybridDecomposition,
        dictionary,
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uspCloud& owner
        ),
        (dict, mesh, owner)
    );


    // Constructors

        //- Construct from components
        uspHybridDecomposition
        (
            const dictionary& dict,
            const polyMesh& mesh,
            uspCloud& owner
        );

    //- Destructor
    virtual ~uspHybridDecomposition() = default;


    //- Selector
    static autoPtr<uspHybridDecomposition> New
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& owner
    );


    // Access

        //- Return the dictionary
        const dictionary& dict() const;


    // Member Functions

        virtual void decompose() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
