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

\*---------------------------------------------------------------------------*/

#include "moleculeIonDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeIonDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    uspReaction,
    moleculeIonDissociationIonisation,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeIonDissociationIonisation::moleculeIonDissociationIonisation
(
    const Time& t,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspReaction(t, cloud, dict),
    dissociationProducts_(),
    ionisationProducts_(),
    heatOfReactionDissociationAB_(),
    heatOfReactionIonisationAB_(),
    nTotABDissociationReactions_(0),
    nTotABIonisationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReationsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeIonDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::moleculeIonDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorInFunction
            << "Reactant molecules cannot be same species."
            << exit(FatalError);
    }

    forAll(reactantIds_, r)
    {
        // check that the first reactant is a 'molecule'

        if (r == 0)
        {
            const auto& constPropReact = cloud_.constProps(reactantIds_[r]);

            label rDof = constPropReact.rotationalDoF();

            if (rDof < 1)
            {
                FatalErrorInFunction
                    << "First reactant must be a molecule (not an atom): "
                    << reactants_[r]
                    << exit(FatalError);
            }

            // check that reactant one only has a single
            // vibrational degree of freedom

            label vDof = constPropReact.vibrationalDoF();

            if (vDof > 1)
            {
                FatalErrorInFunction
                    << "Reactions are currently only implemented "
                    << "for monatomic and diatomic species"
                    << " This is a polyatomic:" << reactants_[r]
                    << exit(FatalError);
            }
        }
    }

    // reading in dissociation products

    wordPair dissociationProducts
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    dissociationProducts_ = labelPair(-1, -1);

    const wordPair& productsForDiss = dissociationProducts;

    forAll(dissociationProducts_, p)
    {
        dissociationProducts_[p] = cloud_.typeIdList().find(productsForDiss[p]);

        if (dissociationProducts_[p] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << productsForDiss[p]
                << exit(FatalError);
        }
    }

    // reading in ionisation products

    wordPair ionisationProducts
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    ionisationProducts_ = labelPair(-1, -1);

    const wordPair& productsForIon = ionisationProducts;

    ionisationProducts_ = labelPair(-1, -1);

    forAll(ionisationProducts_, p)
    {
        ionisationProducts_[p] = cloud_.typeIdList().find(productsForIon[p]);

        if (ionisationProducts_[p] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << productsForIon[p]
                << exit(FatalError);
        }
    }

    // check ionisation product two is an electron
    forAll(ionisationProducts_, p)
    {
        if (p == 1)
        {
            const label charge =
                cloud_.constProps(ionisationProducts_[p]).charge();

            if (charge != -1)
            {
                FatalErrorInFunction
                    << "Second ionisation product must be an electron: "
                    << productsForIon[p]
                    << exit(FatalError);
            }
        }
    }

    heatOfReactionDissociationAB_ =
        propsDict_.get<scalar>("heatOfReactionDissociation");
    heatOfReactionIonisationAB_ =
        propsDict_.get<scalar>("heatOfReactionIonisation");
}


bool Foam::moleculeIonDissociationIonisation::tryReactMolecules
(
    label typeIdP,
    const label typeIdQ
) const
{
    const label reactantPId = reactantIds_.find(typeIdP);
    const label reactantQId = reactantIds_.find(typeIdQ);

    if
    (
        (reactantPId != reactantQId)
     && (reactantPId != -1)
     && (reactantQId != -1)
    )
    {
        return true;
    }

    return false;
}


void Foam::moleculeIonDissociationIonisation::reaction
(
    uspParcel& p,
    uspParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    if (typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        react
        (
            p,
            q          
        );
    }
    else if (p.typeId() == reactantIds_[1] && q.typeId() == reactantIds_[0])
    {
        react
        (
            q,
            p          
        );
    }
}

void Foam::moleculeIonDissociationIonisation::react
(
    uspParcel& p,
    uspParcel& q
)
{
    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(2, 0.0);

    relax_ = true;

    bool dissocReactionP = false;
    bool ionisationReactionP = false;

    // 2 reactions possible

    // 1. Dissociation of P
    // 2. Ionisation of P

    scalar EcPP = 0.0;
    label idP = cloud_.constProps(p.typeId()).charDissQuantumLevel()[0];
    label imaxP = 0;

    // calculate if a dissociation of species P is possible
    EcPP = translationalEnergy(p, q) + EVib(p);

    imaxP = EcPP/(physicoChemical::k.value()*thetaV(p));

    if (imaxP - idP > 0)
    {
        // Dissociation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[0] = 1.0;
    }

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of species P is possible
    EcPP = translationalEnergy(p, q) + EEle(p);

    if ((EcPP - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm(reactionProbabilities.size(),
                                                        0.0);
        scalar sumProbability = 0.0;

        probNorm = reactionProbabilities
                /totalReactionProbability;

        forAll(probNorm, i)
        {
            // If current reaction can't occur, don't check for it
            if (probNorm[i] > VSMALL)
            {
                sumProbability += probNorm[i];

                if (sumProbability > cloud_.rndGen().sample01<scalar>())
                {
                    // Current reaction is to occur

                    if (i == 0)
                    {
                        // Dissociation is to occur
                        dissocReactionP = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Ionisation reaction is to occur
                        ionisationReactionP = true;
                        break;
                    }
                }
            }
        }
    }

    if (dissocReactionP)
    {
        ++nDissociationReactionsPerTimeStep_;
        ++nTotABDissociationReactions_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDissociationAB_,
                dissociationProducts_,
                p,
                q
            );
        }
    }

    if (ionisationReactionP)
    {
        ++nIonisationReationsPerTimeStep_;
        ++nTotABIonisationReactions_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonisationAB_,
                ionisationProducts_,
                p,
                q
            );
        }
    }
}


bool Foam::moleculeIonDissociationIonisation::outputResults
(
    const label counterIndex
)
{
    bool write = uspReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA =
        cloud_.typeIdList()[dissociationProducts_[0]];
    const word& productB =
        cloud_.typeIdList()[dissociationProducts_[1]];

    const word& productC =
        cloud_.typeIdList()[ionisationProducts_[0]];
    const word& productD =
        cloud_.typeIdList()[ionisationProducts_[1]];

    if (write)
    {
        // measure density
        const auto& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        List<label> mols (2, 0);

        forAll(cellOccupancy, c)
        {
            const auto& parcelsInCell = cellOccupancy[c];

            for (uspParcel* p : parcelsInCell)
            {
                label id = reactantIds_.find(p->typeId());

                if (id != -1)
                {
                    ++mols[id];
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;

            reactionRate1 =
                (nTotABDissociationReactions* cloud_.nParticle())
               /(counterIndex*deltaT*numberDensities_[0]*numberDensities_[1]*volume);

            reactionRate2 =
                (nTotABIonisationReactions*cloud_.nParticle())
               /(counterIndex*deltaT*numberDensities_[0]*numberDensities_[1]*volume);

            Info<< "Dissociation type I reaction "
                <<  reactantA
                << " + " << reactantB << " --> "
                << productA << " + "
                << productB << " + "
                << reactantB << ", reaction rate = "
                << reactionRate1  << nl
                << "Ionisation type I reaction "
                <<  reactantA << " + "
                << reactantB << " --> "
                << productC << " + "
                << productD << " + "
                << reactantB << ", reaction rate = "
                << reactionRate2
                << endl;
        }
    }
    else
    {
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;

        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReationsPerTimeStep =
            nIonisationReationsPerTimeStep_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReationsPerTimeStep, sumOp<label>());
        }

        if (nTotABDissociationReactions > VSMALL)
        {
            Info << "Dissociation type I reaction " <<  reactantA
                << " + " << reactantB << " --> "
                << productA << " + "
                << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotABIonisationReactions > VSMALL)
        {
            Info << "Ionisation type I reaction " <<  reactantA
                << " + " << reactantB << " --> "
                << productC << " + "
                << productD << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReationsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReationsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
