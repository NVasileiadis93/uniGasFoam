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

#include "moleculeElectronDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeElectronDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    uspReaction,
    moleculeElectronDissociationIonisation,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeElectronDissociationIonisation::
moleculeElectronDissociationIonisation
(
    const Time& t,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspReaction(t, cloud, dict),
    productIdsDiss_(),
    productIdsIon_(),
    heatOfReactionDiss_(),
    heatOfReactionIon_(),
    nDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeElectronDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::moleculeElectronDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorInFunction
            << "Reactant molecules cannot be same species."
            << exit(FatalError);
    }

    // check that reactant one is a 'molecule''

    label rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDoF();

    if (rDof1 < 1)
    {
        FatalErrorInFunction
            << "First reactant must be a molecule "
            << "(not an atom or an electron): " << reactants_[0]
            << exit(FatalError);
    }

    label vDof = cloud_.constProps(reactantIds_[0]).vibrationalDoF();

    if (vDof > 1)
    {
        FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    // check that reactant two is an 'electron'

    label charge = cloud_.constProps(reactantIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second reactant must be an electron, not "
            << reactants_[1]
            << exit(FatalError);
    }

    // reading in dissociation products

    wordPair productMoleculesDissociation
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    productIdsDiss_ = labelPair(-1, -1);

    forAll(productMoleculesDissociation, r)
    {
        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] =
                cloud_.typeIdList().find(productMoleculesDissociation[r]);

            // check that reactants belong to the typeIdList
            if (productIdsDiss_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << productMoleculesDissociation[r]
                    << exit(FatalError);
            }
        }

        // check that product one is an 'atom'

        label rDof3 = cloud_.constProps(productIdsDiss_[0]).rotationalDoF();

        if (rDof3 != 0)
        {
            FatalErrorInFunction
                << "First product must be an atom (not a molecule): "
                << productMoleculesDissociation[0]
                << exit(FatalError);
        }

        // check that product two is an 'atom'

        label rDof4 = cloud_.constProps(productIdsDiss_[1]).rotationalDoF();

        if (rDof4 != 0)
        {
            FatalErrorInFunction
                << "Second product must be an atom (not a molecule): "
                << productMoleculesDissociation[1]
                << exit(FatalError);
        }
    }

    // reading in ionisation products

    wordPair productMoleculesIonisation
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    productIdsIon_ = labelPair(-1 , -1);

    forAll(productMoleculesIonisation, r)
    {
        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] =
                cloud_.typeIdList().find(productMoleculesIonisation[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << productMoleculesIonisation[r]
                    << exit(FatalError);
            }
        }

        // check that product one is a 'molecule'

        label rDof5 = cloud_.constProps(productIdsIon_[0]).rotationalDoF();

        if (rDof5 < 1)
        {
            FatalErrorInFunction
                << "First product must be a molecule ("
                << "not an atom): "
                << productMoleculesIonisation[0]
                << exit(FatalError);
        }

        // check that product two is an 'electron'

        const label charge = cloud_.constProps(productIdsIon_[1]).charge();

        if (charge != -1)
        {
            FatalErrorInFunction
                << "Second product must be an electron: "
                << productMoleculesIonisation[1]
                << exit(FatalError);
        }
    }

    heatOfReactionDiss_ = propsDict_.get<scalar>("heatOfReactionDissociation");
    heatOfReactionIon_ = propsDict_.get<scalar>("heatOfReactionIonisation");
}


bool Foam::moleculeElectronDissociationIonisation::tryReactMolecules
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


void Foam::moleculeElectronDissociationIonisation::reaction
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

void Foam::moleculeElectronDissociationIonisation::react
(
    uspParcel& p,
    uspParcel& q
)
{
    relax_ = true;

    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(2, 0.0);

    bool dissocReaction = false;
    bool ionisationReaction = false;

    scalar EcPPIon = 0.0;
    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of P is possible
    EcPPIon = translationalEnergy(p, q) + EEle(p);

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        totalReactionProbability += 1.0;
        reactionProbabilities[0] = 1.0;
    }

    scalar EcPPDiss = 0.0;
    label idP = charDissLevel(p);
    label imaxP = 0;

    // calculate if a dissociation of P is possible
    EcPPDiss = translationalEnergy(p, q) + EVib(p);

    imaxP = EcPPDiss/(physicoChemical::k.value()*thetaV(p));

    if (imaxP - idP > 0)
    {
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

        forAll(probNorm, i)
        {
            probNorm[i] = reactionProbabilities[i]
                                            /totalReactionProbability;

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
                        ionisationReaction = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Dissociation reaction is to occur
                        dissocReaction = true;
                        break;
                    }
                }
            }
        }
    }

    // Perform a dissociation reaction
    if (dissocReaction)
    {
        nDissociationReactions_++;
        nDissociationReactionsPerTimeStep_++;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDiss_,
                productIdsDiss_,
                p,
                q
            );
        }
    }

    // Perform an ionisation reaction

    if (ionisationReaction)
    {
        ++nIonisationReactions_;
        ++nIonisationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIon_,
                productIdsIon_,
                p,
                q
            );
        }
    }
}


bool Foam::moleculeElectronDissociationIonisation::outputResults
(
    const label counterIndex
)
{
    bool write = uspReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIdsDiss_[0]];
    const word& productB = cloud_.typeIdList()[productIdsDiss_[1]];

    const word& productC = cloud_.typeIdList()[productIdsIon_[0]];
    const word& productD = cloud_.typeIdList()[productIdsIon_[1]];

    if (write)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        List<label> mols (2, 0);
        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            const auto& parcelsInCell = cellOccupancy[c];

            for (uspParcel* p : parcelsInCell)
            {
                const label id = reactantIds_.find(p->typeId());

                if (id != -1)
                {
                    ++mols[id];
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotReactionsDissociation = nDissociationReactions_;
        label nTotReactionsIonisation = nIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateDissociation = 0.0;

            reactionRateDissociation =
                (nTotReactionsDissociation*cloud_.nParticle())
               /(
                    counterIndex
                   *deltaT
                   *numberDensities_[0]
                   *numberDensities_[1]
                   *volume
                );

            Info<< "Electron dissociation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + " << reactantB
                << ", reaction rate = " << reactionRateDissociation
                << endl;

            scalar reactionRateIonisation = 0.0;

            reactionRateIonisation =
                (nTotReactionsIonisation*cloud_.nParticle())
               /(
                    counterIndex
                   *deltaT
                   *numberDensities_[0]
                   *numberDensities_[1]
                   *volume
                );

            Info<< "Electron ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productC << " + " << productD << " + " << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotReactionsDissociation = nDissociationReactions_;
        label nTotReactionsIonisation = nIonisationReactions_;
        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());

            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotReactionsDissociation > VSMALL)
        {
            Info<< "Electron dissociation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIonisation > VSMALL)
        {
            Info<< "Electron ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
