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
    forwardAssociativeIonisationDissimilarSpecies

Description


SourceFiles
    forwardAssociativeIonisationDissimilarSpecies.C


\*---------------------------------------------------------------------------*/

#ifndef forwardAssociativeIonisationDissimilarSpecies_H
#define forwardAssociativeIonisationDissimilarSpecies_H

#include "uniGasReaction.H"
#include "vector.H"
#include "uniGasCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class forwardAssociativeIonisationDissimilarSpecies Declaration
\*---------------------------------------------------------------------------*/

class forwardAssociativeIonisationDissimilarSpecies
:
    public uniGasReaction
{
    // Private data

        label intermediateId_;

        labelPair associativeIonisationProductIds_;
        labelPair ionisationPProductIds_;
        labelPair ionisationQProductIds_;

        scalar heatOfReactionRecombination_;
        scalar heatOfReactionIntermediateIonisation_;
        scalar heatOfReactionIonisationP_;
        scalar heatOfReactionIonisationQ_;

        label nTotalIonisationPReactions_;
        label nTotalIonisationQReactions_;
        label nTotalAssociativeIonisationReactions_;

        label nIonisationPReactionsPerTimeStep_;
        label nIonisationQReactionsPerTimeStep_;
        label nAssociativeIonisationReactionsPerTimeStep_;


    // Private Member Functions

        void setProperties();
        
        void react
        (
            uniGasParcel& p,
            uniGasParcel& q
        );


public:

    //- Runtime type information
    TypeName("forwardAssociativeIonisationDissimilarSpecies");

    // Constructors

        //- Construct from components
        forwardAssociativeIonisationDissimilarSpecies
        (
            const Time& t,
            uniGasCloud& cloud,
            const dictionary& dict
        );

        //- Destructor
        virtual ~forwardAssociativeIonisationDissimilarSpecies() = default;

    // Member Functions

        //- Initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- Apply a reaction between parcel p and q
        void reaction(uniGasParcel& p, uniGasParcel& q);

        //- Return bool if molecules belong to this reaction model
        bool tryReactMolecules(label typeIdP, const label typeIdQ) const;

        bool outputResults(const label counterIndex);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
