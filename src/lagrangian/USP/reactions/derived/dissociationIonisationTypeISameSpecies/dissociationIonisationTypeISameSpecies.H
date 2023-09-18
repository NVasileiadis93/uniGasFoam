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
    dissociationIonisationTypeISameSpecies

Description


SourceFiles
    dissociationIonisationTypeISameSpecies.C


\*---------------------------------------------------------------------------*/

#ifndef dissociationIonisationTypeISameSpecies_H
#define dissociationIonisationTypeISameSpecies_H

#include "uspReaction.H"
#include "vector.H"
#include "uspCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class dissociationIonisationTypeISameSpecies Declaration
\*---------------------------------------------------------------------------*/

class dissociationIonisationTypeISameSpecies
:
   public uspReaction
{
    // Private data

        labelPair productIdsIonisation_;
        labelPair productIdsDissociation_;

        label nTotIonisationReactions_;
        label nIonisationReactionsPerTimeStep_;
        label nTotDissociationReactions_;
        label nDissociationReactionsPerTimeStep_;

        scalar heatOfReactionIonisation_;
        scalar heatOfReactionDissociation_;


    // Private Member Functions

        void setProperties();
        
        void react
        (
            uspParcel& p,
            uspParcel& q
        );


public:

    //- Runtime type information
    TypeName("dissociationIonisationTypeISameSpecies");

    // Constructors

        //- Construct from components
        dissociationIonisationTypeISameSpecies
        (
            const Time& t,
            uspCloud& cloud,
            const dictionary& dict
        );

        //- Destructor
        virtual ~dissociationIonisationTypeISameSpecies() = default;


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
