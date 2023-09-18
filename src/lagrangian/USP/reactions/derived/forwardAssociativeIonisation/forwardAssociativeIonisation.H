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
    forwardAssociativeIonisation

Description


SourceFiles
    forwardAssociativeIonisation.C


\*---------------------------------------------------------------------------*/

#ifndef forwardAssociativeIonisation_H
#define forwardAssociativeIonisation_H

#include "uspReaction.H"
#include "vector.H"
#include "uspCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class forwardAssociativeIonisation Declaration
\*---------------------------------------------------------------------------*/

class forwardAssociativeIonisation
:
   public uspReaction
{
    // Private data

        label intermediateId_;

        labelPair associativeIonisationProductIds_;
        labelPair ionisationProductIds_;

        scalar heatOfReactionRecombination_;
        scalar heatOfReactionIntermediateIonisation_;
        scalar heatOfReactionIonisation_;

        label nTotalIonisationReactions_;
        label nTotalAssociativeIonisationReactions_;

        label nIonisationReactionsPerTimeStep_;
        label nAssociativeIonisationReactionsPerTimeStep_;


    // Private Member Functions

        void setProperties();
        
        void react
        (
            uspParcel& p,
            uspParcel& q
        );


public:

    //- Runtime type information
    TypeName("forwardAssociativeIonisation");

    // Constructors

        //- Construct from components
        forwardAssociativeIonisation
        (
            const Time& t,
            uspCloud& cloud,
            const dictionary& dict
        );

        //- Destructor
        virtual ~forwardAssociativeIonisation() = default;


    // Member Functions

        //- Initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- Apply a reaction between parcel p and q
        void reaction(uspParcel& p, uspParcel& q);

        //- Return true if molecules belong to this reaction model
        bool tryReactMolecules(label typeIdP, const label typeIdQ) const;

        bool outputResults(const label counterIndex);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
