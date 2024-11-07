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
    Foam::uniGasReactions

Description
    Stores all the information for the uniGasReactions models defined within
    the uniGasReactionsDict, and selects & builds the models automatically.

Source Files
    uniGasReactions.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasReactions_H
#define uniGasReactions_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "uniGasReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class uniGasReactions Declaration
\*---------------------------------------------------------------------------*/

class uniGasReactions
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary chemReactDict_;

        //- Number of uniGasReactions
        label nReactions_;

        PtrList<entry> reactionsList_;
        List<word> reactionNames_;
        List<label> reactionIds_;
        List<autoPtr<uniGasReaction>> reactions_;

        //- Given two pair typeIds, one is returned the modelId number
        //  in the reactions_ data member
        List<List<label>> pairAddressing_;

        label counter_;


public:

    // Constructors

        //- Null Constructor
        uniGasReactions
        (
            const Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        uniGasReactions
        (
            const Time& t,
            const polyMesh& mesh,
            uniGasCloud& cloud
        );


        virtual ~uniGasReactions() = default;


    // Member Functions

        //- Initial configuration
        void initialConfiguration();

        void outputData();


        // Access

            inline label returnModelId(uniGasParcel& p, uniGasParcel& q);

            inline List<autoPtr<uniGasReaction>>& reactions();
            inline const List<autoPtr<uniGasReaction>>& reactions() const;

            inline const List<List<label>>& pairModelAddressing() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniGasReactionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
