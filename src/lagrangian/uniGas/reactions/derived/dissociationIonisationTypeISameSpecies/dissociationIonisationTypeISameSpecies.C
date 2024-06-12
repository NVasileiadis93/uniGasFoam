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

\*---------------------------------------------------------------------------*/

#include "dissociationIonisationTypeISameSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationTypeISameSpecies, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    dissociationIonisationTypeISameSpecies,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissociationIonisationTypeISameSpecies::
dissociationIonisationTypeISameSpecies
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    productIdsIonisation_(),
    productIdsDissociation_(),
    nTotIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    nTotDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    heatOfReactionIonisation_(),
    heatOfReactionDissociation_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissociationIonisationTypeISameSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::dissociationIonisationTypeISameSpecies::setProperties()
{

    if (reactantIds_[0] != reactantIds_[1])
    {
        FatalErrorInFunction
            << "Reactants must be same species!"
            << exit(FatalError);
    }

    if (rDof1_ < VSMALL)
    {
        FatalErrorInFunction
            << "Reactants must be molecules: "
            << reactants_[0]
            << exit(FatalError);
    }

    if (vDof1_ > 1)
    {
         FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species."
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    // reading in ionisation products

    const wordPair productMoleculesIonisation
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    productIdsIonisation_ = labelPair(-1, -1);

    forAll(productIdsIonisation_, r)
    {
        productIdsIonisation_[r] =
            cloud_.typeIdList().find(productMoleculesIonisation[r]);

        // check that products belong to the typeIdList
        if (productIdsIonisation_[r] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << productMoleculesIonisation[r]
                << exit(FatalError);
        }
    }

    label rDof2 = cloud_.constProps(productIdsIonisation_[0]).rotationalDoF();

    if (rDof2 < VSMALL)
    {
        FatalErrorInFunction
            << "First ionisation product must be a molecular ion: "
            << productMoleculesIonisation[0]
            << exit(FatalError);
    }

    const label charge = cloud_.constProps(productIdsIonisation_[1]).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second ionisation product must be an electron: "
            << productMoleculesIonisation[1]
            << exit(FatalError);
    }

    // reading in dissociation products

    const wordPair productMoleculesDissociation
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    productIdsDissociation_ = labelPair(-1, -1);

    forAll(productIdsDissociation_, r)
    {
        productIdsDissociation_[r] =
            cloud_.typeIdList().find(productMoleculesDissociation[r]);

        // check that products belong to the typeIdList
        if (productIdsDissociation_[r] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << productMoleculesDissociation[r]
                << exit(FatalError);
        }

        // check that products are 'atoms'

        label iD = productIdsDissociation_[r];

        label rDof3 = cloud_.constProps(iD).rotationalDoF();

        if (rDof3 > 1)
        {
            FatalErrorInFunction
                << "Dissociation product must be an atom (not a molecule): "
                << productMoleculesDissociation[r]
                << exit(FatalError);
        }
    }

    heatOfReactionIonisation_ =
        propsDict_.get<scalar>("heatOfReactionIonisation");

    heatOfReactionDissociation_ =
        propsDict_.get<scalar>("heatOfReactionDissociation");
}


bool Foam::dissociationIonisationTypeISameSpecies::tryReactMolecules
(
    label typeIdP,
    const label typeIdQ
) const
{
    const label reactantPId = reactantIds_.find(typeIdP);
    const label reactantQId = reactantIds_.find(typeIdQ);

    if
    (
        (reactantPId == reactantQId)
     && (reactantPId != -1)
     && (reactantQId != -1)
    )
    {
        return true;
    }

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


void Foam::dissociationIonisationTypeISameSpecies::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    if (typeIdP == typeIdQ && typeIdP == reactantIds_[0])
    {
        react
        (
            p,
            q          
        );
    }
}

void Foam::dissociationIonisationTypeISameSpecies::react
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

    EcPP = translationalEnergy(p, q) + EVib(p);
    label imaxP = EcPP/(physicoChemical::k.value()*thetaV(p));
    label idP = thetaD(p)/thetaV(p);

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
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }

    scalar EcPQ = translationalEnergy(p, q) + EVib(q);

    label imaxQ = EcPQ/(physicoChemical::k.value()*thetaV(q));
    label idQ = thetaD(q)/thetaV(q);

    if (imaxQ - idQ > 0)
    {
        // Dissociation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[2] = 1.0;
    }

    scalar ionisationEnergyQ =
        cloud_.constProps(q.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of Q is possible
    EcPQ = translationalEnergy(p, q) + EEle(q);

    if ((EcPQ - ionisationEnergyQ) > VSMALL)
    {
        totalReactionProbability += 1.0;
        reactionProbabilities[3] = 1.0;
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm(reactionProbabilities.size(),                                                                   0.0);
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
                        // Ionisation is to occur
                        dissocReactionP = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Dissociation reaction is to occur
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

    // dissociation of P
    if (dissocReactionP)
    {
        ++nTotDissociationReactions_;
        ++nDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDissociation_,
                productIdsDissociation_,
                p,
                q
            );
        }
    }

    // dissociation of Q
    if (dissocReactionQ )
    {
        ++nTotDissociationReactions_;
        ++nDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateQ
            (
                heatOfReactionDissociation_,
                productIdsDissociation_,
                p,
                q
            );
        }
    }

    // ionisation of P
    if (ionisationReactionP)
    {
        ++nTotIonisationReactions_;
        ++nIonisationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonisation_,
                productIdsIonisation_,
                p,
                q
            );
        }
    }

    // ionisation of Q
    if (ionisationReactionQ)
    {
        ++nTotIonisationReactions_;
        ++nIonisationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseQ
            (
                heatOfReactionIonisation_,
                productIdsIonisation_,
                p,
                q
            );
        }
    }
}


bool Foam::dissociationIonisationTypeISameSpecies::outputResults
(
    const label counterIndex
)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIdsDissociation_[0]];
    const word& productB = cloud_.typeIdList()[productIdsDissociation_[1]];

    const word& productC = cloud_.typeIdList()[productIdsIonisation_[0]];
    const word& productD = cloud_.typeIdList()[productIdsIonisation_[1]];

    if (write)
    {
        // measure density

        const List<DynamicList<uniGasParcel*>>& cellOccupancy
            = cloud_.cellOccupancy();

        volume_ = 0.0;

        label molsReactants = 0;

        forAll(cellOccupancy, c)
        {
            const List<uniGasParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                uniGasParcel* p = parcelsInCell[pIC];

                if (reactantIds_.find(p->typeId()) != -1)
                {
                    ++molsReactants;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotReactions = nTotReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(molsReactants, sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotReactions, sumOp<label>());
        }

        numberDensities_[0] = (molsReactants*cloud().nParticle())/volume;
        numberDensities_[1] = (molsReactants*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateDissociation = 0.0;

            reactionRateIonisation =
                (
                    nTotIonisationReactions_
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Ionisation type I reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productC << " + " << productD << " + " << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            reactionRateDissociation =
                (
                    nTotDissociationReactions_
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Dissociation type I reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + " << reactantB
                << ", reaction rate = " << reactionRateDissociation
                << endl;
        }
    }
    else
    {
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactions = nTotIonisationReactions_;
        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotDissociationReactions > VSMALL)
        {
            Info<< "Dissociation type I reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + " << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactions > VSMALL)
        {
            Info<< "Ionisation type I reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + " << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0;
    nIonisationReactionsPerTimeStep_ = 0;

    return write;
}


// ************************************************************************* //
