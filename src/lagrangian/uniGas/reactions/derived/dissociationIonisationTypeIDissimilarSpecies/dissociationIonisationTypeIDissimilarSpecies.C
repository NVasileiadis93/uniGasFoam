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

\*---------------------------------------------------------------------------*/

#include "dissociationIonisationTypeIDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationTypeIDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    dissociationIonisationTypeIDissimilarSpecies,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissociationIonisationTypeIDissimilarSpecies::
dissociationIonisationTypeIDissimilarSpecies
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    dissociationProducts_(),
    ionisationProducts_(),
    heatOfReactionDissociationAB_(),
    heatOfReactionIonisationAB_(),
    heatOfReactionDissociationCD_(),
    heatOfReactionIonisationCD_(),
    nTotABDissociationReactions_(0),
    nABDissociationReactionsPerTimeStep_(0),
    nTotABIonisationReactions_(0),
    nABIonisationReactionsPerTimeStep_(0),
    nTotCDDissociationReactions_(0),
    nCDDissociationReactionsPerTimeStep_(0),
    nTotCDIonisationReactions_(0),
    nCDIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissociationIonisationTypeIDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::dissociationIonisationTypeIDissimilarSpecies::setProperties()
{
    // check that reactants are 'molecules'

    if (rDof1_ < 1)
    {
        FatalErrorInFunction
            << "Reactant must be a molecule (not an atom): "
            << reactants_[0]
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
        FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    if (rDof2_ < 1)
    {
        FatalErrorInFunction
            << "Reactant must be a molecule (not an atom): "
            << reactants_[1]
            << exit(FatalError);
    }

    if (vDof2_ > 1)
    {
        FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[1]
            << exit(FatalError);
    }

    // reading in dissociation products

    List<wordPair> dissociationProducts
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    if (dissociationProducts.size() != reactantIds_.size())
    {
        FatalErrorInFunction
            << "number of reactant molecules to be dissociated = "
            << reactantIds_.size()
            << " is not the same as the number of products = "
            << dissociationProducts.size()
            << exit(FatalError);
    }


    dissociationProducts_.setSize(dissociationProducts.size());

    forAll(dissociationProducts_, r)
    {
        const wordPair& productsForDiss = dissociationProducts[r];

        dissociationProducts_[r] = labelPair(-1, -1);

        forAll(dissociationProducts_[r], p)
        {
            dissociationProducts_[r][p] =
                cloud_.typeIdList().find(productsForDiss[p]);

            if (dissociationProducts_[r][p] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: " << productsForDiss[p] << nl
                    << exit(FatalError);
            }

            // check that products are 'atoms'
            label iD = dissociationProducts_[r][p];

            label rDof3 = cloud_.constProps(iD).rotationalDoF();

            if (rDof3 > 1)
            {
                FatalErrorInFunction
                    << "Dissociation product must be an atom (not a molecule): "
                    << productsForDiss[r][p]
                    << nl
                    << exit(FatalError);
            }
        }
    }

    // reading in ionisation products

    List<wordPair> ionisationProducts
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    if (ionisationProducts.size() != reactantIds_.size())
    {
        FatalErrorInFunction
            << "number of reactant molecules to be ionised = "
            << reactantIds_.size()
            << " is not the same as the number of products = "
            << ionisationProducts.size()
            << exit(FatalError);
    }


    ionisationProducts_.setSize(ionisationProducts.size());

    forAll(ionisationProducts_, r)
    {
        const wordPair& productsForIon = ionisationProducts[r];

        ionisationProducts_[r] = labelPair(-1, -1);

        forAll(ionisationProducts_[r], p)
        {
            ionisationProducts_[r][p] =
                cloud_.typeIdList().find(productsForIon[p]);

            if (ionisationProducts_[r][p] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: " << productsForIon[p]
                    << exit(FatalError);
            }
        }

        // check ionisation product two is an 'electron'
        forAll(ionisationProducts_[r], p)
        {
            if (p == 1)
            {
                const label charge =
                    cloud_.constProps(ionisationProducts_[r][p]).charge();

                if (charge != -1)
                {
                    FatalErrorInFunction
                        << "Second ionisation product must be an electron: "
                        << productsForIon[p]
                        << exit(FatalError);
                }
            }
        }
    }

    heatOfReactionDissociationAB_ =
        propsDict_.get<scalar>("heatOfReactionDissociationAB");
    heatOfReactionIonisationAB_ =
        propsDict_.get<scalar>("heatOfReactionIonisationAB");
    heatOfReactionDissociationCD_ =
        propsDict_.get<scalar>("heatOfReactionDissociationCD");
    heatOfReactionIonisationCD_ =
        propsDict_.get<scalar>("heatOfReactionIonisationCD");
}


bool Foam::dissociationIonisationTypeIDissimilarSpecies::tryReactMolecules
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


void Foam::dissociationIonisationTypeIDissimilarSpecies::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
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

void Foam::dissociationIonisationTypeIDissimilarSpecies::react
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    relax_ = true;
    
    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(4, 0.0);

    bool dissocReactionP = false;
    bool ionisationReactionP = false;
    bool dissocReactionQ = false;
    bool ionisationReactionQ = false;

    // 4 reactions possible

    // 1. Dissociation of P
    // 2. Ionisation of P
    // 3. Dissociation of Q
    // 4. Ionisation of Q

    scalar EcPP = 0.0;
    label idP = charDissLevel(p);
    label imaxP = 0;

    // calculate if a dissociation of P is possible
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

    // calculate if an ionisation of P is possible
    EcPP = translationalEnergy(p, q) + EEle(p);

    if ((EcPP - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }

    scalar EcQP = 0.0;
    label idQ = charDissLevel(q);
    label imaxQ = 0;

    // calculate if a dissociation of Q is possible
    EcQP = translationalEnergy(p, q) + EVib(q);

    imaxQ = EcQP/(physicoChemical::k.value()*thetaV(q));

    if (imaxQ - idQ > 0)
    {
        // Dissociation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[2] = 1.0;
    }

    ionisationEnergy =
        cloud_.constProps(q.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of Q is possible
    EcQP = translationalEnergy(p, q) + EEle(q);

    if ((EcQP - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[3] = 1.0;
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm(reactionProbabilities.size(),
                                                                0.0);
        scalar sumProbability = 0.0;

        probNorm =
            reactionProbabilities/totalReactionProbability;

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
                    if (i == 2)
                    {
                        // Dissociation is to occur
                        dissocReactionQ = true;
                        break;
                    }
                    if (i == 3)
                    {
                        // Ionisation reaction is to occur
                        ionisationReactionQ = true;
                        break;
                    }
                }
            }
        }
    }

    if (dissocReactionP)
    {
        ++nTotABDissociationReactions_;
        ++nABDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDissociationAB_,
                dissociationProducts_[0],
                p,
                q
            );
        }
    }

    if (dissocReactionQ)
    {
        ++nTotCDDissociationReactions_;
        ++nCDDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateQ
            (
                heatOfReactionDissociationCD_,
                dissociationProducts_[1],
                p,
                q
            );
        }
    }

    if (ionisationReactionP)
    {
        nTotABIonisationReactions_++;
        nABIonisationReactionsPerTimeStep_++;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonisationAB_,
                ionisationProducts_[0],
                p,
                q
            );
        }
    }

    if (ionisationReactionQ)
    {
        ++nTotCDIonisationReactions_;
        ++nCDIonisationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseQ
            (
                heatOfReactionIonisationCD_,
                ionisationProducts_[1],
                p,
                q
            );
        }
    }
}


bool Foam::dissociationIonisationTypeIDissimilarSpecies::outputResults
(
    const label counterIndex
)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA =
        cloud_.typeIdList()[dissociationProducts_[0][0]];
    const word& productB =
        cloud_.typeIdList()[dissociationProducts_[0][1]];
        
    const word& productC =
        cloud_.typeIdList()[dissociationProducts_[1][0]];
    const word& productD =
        cloud_.typeIdList()[dissociationProducts_[1][1]];

    const word& productE =
        cloud_.typeIdList()[ionisationProducts_[0][0]];
    const word& productF =
        cloud_.typeIdList()[ionisationProducts_[0][1]];
        
    const word& productG =
        cloud_.typeIdList()[ionisationProducts_[1][0]];
    const word& productH =
        cloud_.typeIdList()[ionisationProducts_[1][1]];

    if (write)
    {
        // measure density
        const auto& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        List<label> mols (2, 0);

        forAll(cellOccupancy, c)
        {
            const List<uniGasParcel*>& parcelsInCell = cellOccupancy[c];

            for (uniGasParcel* p : parcelsInCell)
            {
                label id = reactantIds_.find(p->typeId());

                if (id != -1)
                {
                    mols[id]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;

            reactionRate1 =
                (
                    nTotABDissociationReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            reactionRate2 =
                (
                    nTotCDDissociationReactions
                   *cloud_.nParticle()
                )
               /(
                   counterIndex*deltaT*numberDensities_[0]
                  *numberDensities_[1]*volume
                );

            reactionRate3 =
                (
                    nTotABIonisationReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            reactionRate4 =
            (
                nTotCDIonisationReactions
               *cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]
            * numberDensities_[1]*volume);

            Info << "Dissociation type I reaction " <<  reactantA
                 << " + " << reactantB << " --> "
                 << productA << " + "
                 << productB << " + " << reactantB
                 << ", reaction rate = " << reactionRate1
                 << nl
                 << "Dissociation type I reaction " <<  reactantB
                 << " + " << reactantA << " --> "
                 << productC << " + "
                 << productD << " + " << reactantA
                 << ", reaction rate = " << reactionRate2
                 << nl
                 << "Ionisation type I reaction " <<  reactantA
                 << " + " << reactantB << " --> "
                 << productE << " + "
                 << productF << " + " << reactantB
                 << ", reaction rate = " << reactionRate3
                 << nl
                 << "Ionisation type I reaction " <<  reactantB
                 << " + " << reactantA << " --> "
                 << productG << " + "
                 << productH << " + " << reactantA
                 << ", reaction rate = " << reactionRate4
                 << endl;
        }
    }
    else
    {
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;

        label nABDissociationReactionsPerTimeStep =
            nABDissociationReactionsPerTimeStep_;
        label nCDDissociationReactionsPerTimeStep =
            nCDDissociationReactionsPerTimeStep_;
        label nABIonisationReactionsPerTimeStep =
            nABIonisationReactionsPerTimeStep_;
        label nCDIonisationReactionsPerTimeStep =
            nCDIonisationReactionsPerTimeStep_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());

            reduce(nABDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nABIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotABDissociationReactions > VSMALL)
        {
            Info << "Dissociation type I reaction " <<  reactantA << " + "
                 << reactantB << " --> " <<  productA
                 << " + "  << productA << " + "
                 << reactantB << " is active, nReactions this "
                 << "time step = " << nABDissociationReactionsPerTimeStep
                 << endl;
        }

        if (nTotCDDissociationReactions > VSMALL)
        {
            Info << "Dissociation type I reaction " <<  reactantB << " + "
                 << reactantA << " --> " << productC
                 << " + " << productD << " + "
                 << reactantA << " is active, nReactions this "
                 << "time step = " << nCDDissociationReactionsPerTimeStep
                 << endl;
        }

        if (nTotABIonisationReactions > VSMALL)
        {
            Info << "Ionisation type I reaction " <<  reactantA << " + "
                 << reactantB << " --> " << productE
                 << " + " << productF << " + "
                 << reactantB << " is active, nReactions this "
                 << "time step = " << nABIonisationReactionsPerTimeStep
                 << endl;
        }

        if (nTotCDIonisationReactions > VSMALL)
        {
            Info << "Ionisation type I reaction " <<  reactantB << " + "
                 << reactantA << " --> " << productG
                 << " + " << productH << " + "
                 << reactantA << " is active, nReactions this "
                 << "time step = " << nCDIonisationReactionsPerTimeStep
                 << endl;
        }
    }

    nABDissociationReactionsPerTimeStep_ = 0;
    nCDDissociationReactionsPerTimeStep_ = 0;
    nABIonisationReactionsPerTimeStep_ = 0;
    nCDIonisationReactionsPerTimeStep_ = 0;

    return write;
}


// ************************************************************************* //
