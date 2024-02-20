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
    Foam::uspBoundaries

Description
    Stores all the information for the boundary models.

    Models are defined within the boundariesDict. This class selects and builds
    the models automatically. We assume that boundaries are divided into three:

    1. Patch models:
       - actions on instantaneous collisions with a patch boundary

    2. cyclic boundary models:
       - actions on instantaneous collisions with a cyclic boundary

    3. General boundary models:
       - time-dependent actions on molecules close to the boundary

Source Files
    uspBoundaries.C
    
\*---------------------------------------------------------------------------*/

#ifndef uspBoundaries_H
#define uspBoundaries_H

#include "List.H"

#include "uspCyclicBoundary.H"
#include "uspPatchBoundary.H"
#include "uspGeneralBoundary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class uspBoundaries Declaration
\*---------------------------------------------------------------------------*/

class uspBoundaries
{
    // Private data

        const Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary uspBoundariesDict_;

        //- Number of boundary models
        label nPatchBoundaryModels_;
        label nCyclicBoundaryModels_;
        label nGeneralBoundaryModels_;

        //- Patch models
        PtrList<entry> patchBoundaryList_;
        List<word> patchBoundaryNames_;
        List<label> patchBoundaryIds_;
        List<fileName> pBFixedPathNames_;
        List<autoPtr<uspPatchBoundary>> patchBoundaryModels_;
        List<label> patchToModelId_;

        //- Cyclic boundary models
        PtrList<entry> cyclicBoundaryList_;
        List<word> cyclicBoundaryNames_;
        List<label> cyclicBoundaryIds_;
        List<fileName> cMFixedPathNames_;
        List<autoPtr<uspCyclicBoundary>> cyclicBoundaryModels_;
        List<label> cyclicBoundaryToModelId_;

        //- General boundary models
        PtrList<entry> generalBoundaryList_;
        List<word> generalBoundaryNames_;
        List<label> generalBoundaryIds_;
        List<fileName> gMFixedPathNames_;
        List<autoPtr<uspGeneralBoundary>> generalBoundaryModels_;


public:

    //- Default name of the boundaries dictionary
    static const word dictName;


    // Constructors

        //- Null Constructor
        uspBoundaries(const Time& t, const polyMesh& mesh);

        //- Constructor
        uspBoundaries(const Time& t, const polyMesh& mesh, uspCloud& cloud);

        ~uspBoundaries() = default;

        void setInitialConfig();

        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        void calculateProps();

        void outputResults();

        void checkPatchBoundaryModels(const polyMesh& mesh);
        void checkCyclicBoundaryModels(const polyMesh& mesh);


    // Member Functions

        label nPatchBoundaryModels() const;
        label nCyclicBoundaryModels() const;
        label nGeneralBoundaryModels() const;

        inline List<autoPtr<uspPatchBoundary>>& patchBoundaryModels();
        inline const List<autoPtr<uspPatchBoundary>>&
            patchBoundaryModels() const;
        inline const List<word>& patchBoundaryNames() const;
        inline const List<label>& patchBoundaryIds() const;
        inline const List<label>& patchToModelIds() const;

        inline List<autoPtr<uspCyclicBoundary>>& cyclicBoundaryModels();
        inline const List<autoPtr<uspCyclicBoundary>>&
            cyclicBoundaryModels() const;
        inline const List<word>& cyclicBoundaryNames() const;
        inline const List<label>& cyclicBoundaryIds() const;
        inline const List<label>& cyclicBoundaryToModelIds() const;


        inline List<autoPtr<uspGeneralBoundary>>& generalBoundaryModels();
        inline const List<autoPtr<uspGeneralBoundary>>&
            generalBoundaryModels() const;
        inline const List<word>& generalBoundaryNames() const;
        inline const List<label>& generalBoundaryIds() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
