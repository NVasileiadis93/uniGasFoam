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
    uniGasReaction

Description
    Abstract class for chemical reactions

SourceFiles
    uniGasReactionI.H
    uniGasReaction.C
    uniGasReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasReaction_H
#define uniGasReaction_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "Random.H"
#include "uniGasParcel.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class uniGasReaction Declaration
\*---------------------------------------------------------------------------*/

class uniGasReaction
{
    // Private Member Functions

        //- Dissociate particle P|Q
        void dissociatePQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const;

        //- Ionise particle P|Q
        void ionisePQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the cloud database
        uniGasCloud& cloud_;

        //- Total number of reactions
        label nTotReactions_;

        //- Number of reactions per time step
        label nReactionsPerTimeStep_;

        //- Properties dictionary
        dictionary propsDict_;

        //- Reaction name
        word reactionName_;

        //- Reactant names
        wordPair reactants_;

        //- Reactant IDs
        labelPair reactantIds_;

        label rDof1_;
        label rDof2_;
        label vDof1_;
        label vDof2_;
        label charge1_;
        label charge2_;

        //- Relax flag
        bool relax_;

        // remove if not used for heat bath
        bool allowSplitting_;

        bool writeRatesToTerminal_;

        scalar volume_;

        List<scalar> numberDensities_;


    // Protected Member Functions

        void setCommonReactionProperties();


        //- Dissociate particle P
        void dissociateP
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const
        {
            dissociatePQ(heatOfReaction, productIDs, p, q);
        }

        //- Dissociate particle Q
        void dissociateQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const
        {
            dissociatePQ(heatOfReaction, reverse(productIDs), q, p);
        }

        //- Ionise particle P
        void ioniseP
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const
        {
            ionisePQ(heatOfReaction, productIDs, p, q);
        }

        //- Ionise particle Q
        void ioniseQ
        (
            const scalar& heatOfReaction,
            const labelPair& productIDs,
            uniGasParcel& p,
            uniGasParcel& q
        ) const
        {
            ionisePQ(heatOfReaction, reverse(productIDs), q, p);
        }

        //- Particle vibrational energy
        scalar EVib(const uniGasParcel& p) const;

        //- Particle rotational energy
        scalar ERot(const uniGasParcel& p) const;

        //- Particle electrical energy level
        scalar EEle(const uniGasParcel& p) const;

        //- Particle mass
        scalar m(const uniGasParcel& p) const;

        //- Reduced mass
        scalar mR(const uniGasParcel& p, const uniGasParcel& q) const;

        //- Tranlational energy
        scalar translationalEnergy
        (
            const uniGasParcel& p,
            const uniGasParcel& q
        ) const;

        //- Omega
        scalar omega(const uniGasParcel& p, const uniGasParcel& q) const;

        //- Particle thetaD
        scalar thetaD(const uniGasParcel& p) const;

        //- Particle thetaV
        scalar thetaV(const uniGasParcel& p) const;

        //- Particle Zref
        scalar Zref(const uniGasParcel& p) const;

        //- Particle refTempZv
        scalar refTempZv(const uniGasParcel& p) const;

        //- Species dissociation level
        label charDissLevel(const uniGasParcel& p) const;

        //- Species number of electronic levels
        label jMax(const uniGasParcel& p) const;

        //- Species rotational models
        label rotationalDof(const uniGasParcel& p) const;

        //- Centre of mass velocity for particle pair
        vector Ucm(const uniGasParcel& p, const uniGasParcel& q) const;

        //- Species electronic energies
        scalarList EEList(const uniGasParcel& p) const;

        //- Species degeneracies
        labelList gList(const uniGasParcel& p) const;

        //- Perform associative ionisation
        void associativeIonisation
        (
            const scalar& heatOfReactionIntermediateIonisation,
            const scalar& heatOfReactionRecombination,
            const labelPair& assIonProductIds,
            uniGasParcel& p,
            uniGasParcel& q
        );
        
        //- Perform exchange reaction
        void exchangePQ
        (
            const scalar& heatOfReactionExchJoules,
            const labelPair& exchangeProductIds,
            uniGasParcel& p,
            uniGasParcel& q
        );
        
        //- Perform charge exchange reaction
        void chargeExchangePQ
        (
            const scalar& heatOfReactionExchJoules,
            const labelPair& chargeExchangeProductIds,
            uniGasParcel& p,
            uniGasParcel& q
        );


public:

    //- Runtime type information
    TypeName("uniGasReaction");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            uniGasReaction,
            dictionary,
            (
                const Time& t,
                uniGasCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );

    // Constructors

        //- Construct from components
        uniGasReaction
        (
            const Time& t,
            uniGasCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<uniGasReaction> New
        (
            const Time& t,
            uniGasCloud& cloud,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uniGasReaction() = default;

    // Member Functions

        virtual void initialConfiguration() = 0;

        //- Apply a reaction between parcel p and q
        virtual void reaction(uniGasParcel& p, uniGasParcel& q) = 0;

        //- Return bool for relax
        //  true : run normal collision model
        //  false : apply reaction to particles
        virtual bool relax() const;

        //- Return bool if particles belong to this reaction model
        virtual bool tryReactMolecules
        (
            const label typeIdP,
            const label typeIdQ
        ) const = 0;

        virtual bool outputResults(const label counterIndex);

        //- Return const access to the cloud
        const uniGasCloud& cloud() const;

        //- Return const access to the total number of reactions
        inline label nTotReactions() const;

        //- Return const access to the number of reactions per time step
        inline label nReactionsPerTimeStep() const;

        //- Return access to the number of reactions per time step
        inline label& nReactionsPerTimeStep();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "uniGasReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
