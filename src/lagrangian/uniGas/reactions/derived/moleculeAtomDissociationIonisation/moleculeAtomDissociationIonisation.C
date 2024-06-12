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

#include "moleculeAtomDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(moleculeAtomDissociationIonisation, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    moleculeAtomDissociationIonisation,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeAtomDissociationIonisation::moleculeAtomDissociationIonisation
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    productIdsDiss_(),
    productIdsIon_(),
    productIdsIon2_(),
    heatOfReactionDiss_(),
    heatOfReactionIon_(),
    heatOfReactionIon2_(),
    nTotReactionsDiss_(0),
    nTotReactionsIon_(0),
    nTotReactionsIon2_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPerTimeStep_(0),
    nIonisationReactions2PerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeAtomDissociationIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::moleculeAtomDissociationIonisation::setProperties()
{
    if (reactantIds_[0] == reactantIds_[1])
    {
        FatalErrorInFunction
            << "Reactant molecules cannot be same species."
            << exit(FatalError);
    }

    // check that first reactant is a 'molecule'

    label rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDoF();

    if (rDof1 < 1)
    {
        FatalErrorInFunction
            << "First reactant must be a molecule (not an atom): "
            << reactants_[0]
            << exit(FatalError);
    }

    label vDof1 = cloud_.constProps(reactantIds_[0]).vibrationalDoF();

    if (vDof1 > 1)
    {
        FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactantIds_[0]
            << exit(FatalError);
    }

    // check that second reactant is an 'atom'

    label rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDoF();

    if (rDof2 > VSMALL)
    {
        FatalErrorInFunction
            << "Second reactant must be atom (not a molecule): "
            << reactants_[1]
            << exit(FatalError);
    }

    // reading in products

    wordPair productMoleculesDiss
    (
        propsDict_.lookup("productsOfDissociatedMolecule")
    );

    wordPair productMoleculesIon
    (
        propsDict_.lookup("productsOfIonisedMolecule")
    );

    wordPair productMoleculesIon2
    (
        propsDict_.lookup("productsOfIonisedAtom")
    );

    productIdsDiss_ = labelPair(-1, -1);

    forAll(productMoleculesDiss, r)
    {
        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] =
                cloud_.typeIdList().find(productMoleculesDiss[r]);

            // check that reactants belong to the typeIdList
            if (productIdsDiss_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: " << productMoleculesDiss[r]
                    << exit(FatalError);
            }
        }

        // check that product one is an 'atom'

        label rDof3 = cloud_.constProps(productIdsDiss_[0]).rotationalDoF();

        if (rDof3 != 0)
        {
            FatalErrorInFunction
                << "First dissociation product must be "
                << "an atom (not a molecule): " << productMoleculesDiss[0]
                << exit(FatalError);
        }

        // check that product two is an 'atom'

        label rDof4 = cloud_.constProps(productIdsDiss_[1]).rotationalDoF();

        if (rDof4 != 0)
        {
            FatalErrorInFunction
                << "Second dissociation product must be "
                << "an atom (not a molecule): " << productMoleculesDiss[1]
                << exit(FatalError);
        }
    }


    productIdsIon_ = labelPair(-1, -1);

    forAll(productMoleculesIon, r)
    {
        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] =
                cloud_.typeIdList().find(productMoleculesIon[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: " << productMoleculesIon[r]
                    << exit(FatalError);
            }
        }

        // check that product one is a 'molecule'

        label rDof5 = cloud_.constProps(productIdsIon_[0]).rotationalDoF();

        if (rDof5 < 0)
        {
            FatalErrorInFunction
                << "First ionisation product must be a "
                << "charged molecule (not an atom/electron): "
                << productMoleculesIon[0]
                << exit(FatalError);
        }

        // check that product two is an 'electron'

        const label charge = cloud_.constProps(productIdsIon_[1]).charge();

        if (charge != -1)
        {
            FatalErrorInFunction
                << "Second ionisation product must be an electron: "
                << productMoleculesIon[1]
                << exit(FatalError);
        }
    }


    productIdsIon2_ = labelPair(-1, -1);

    forAll(productMoleculesIon2, r)
    {
        forAll(productIdsIon2_, r)
        {
            productIdsIon2_[r] =
                cloud_.typeIdList().find(productMoleculesIon2[r]);

            // check that reactants belong to the typeIdList
            if (productIdsIon2_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: " << productMoleculesIon2[r]
                    << exit(FatalError);
            }
        }

        // check that product one is an 'atom'

        label rDof6 = cloud_.constProps(productIdsIon2_[0]).rotationalDoF();

        if (rDof6 > 0)
        {
            FatalErrorInFunction
                << "First ionisation product must be a charged atom: "
                << productMoleculesIon[0]
                << exit(FatalError);
        }

        // check that product two is an 'electron'

        const label charge = cloud_.constProps(productIdsIon2_[1]).charge();

        if (charge != -1)
        {
            FatalErrorInFunction
                << "Second ionisation product must be an electron: "
                << productMoleculesIon[1]
                << exit(FatalError);
        }
    }

    heatOfReactionDiss_ =
        propsDict_.get<scalar>("heatOfReactionDissociation");
    heatOfReactionIon_ =
        propsDict_.get<scalar>("heatOfReactionIonisationMolecule");
    heatOfReactionIon2_ =
        propsDict_.get<scalar>("heatOfReactionIonisationAtom");
}


bool Foam::moleculeAtomDissociationIonisation::tryReactMolecules
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


void Foam::moleculeAtomDissociationIonisation::reaction
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

void Foam::moleculeAtomDissociationIonisation::react
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    relax_ = true;
    
    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(3, 0.0);

    bool dissocReaction = false;
    bool ionisationReaction = false;
    bool atomIonisationReaction = false;

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

    ionisationEnergy =
        cloud_.constProps(q.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of Q is possible
    EcPPIon = translationalEnergy(p, q) + EEle(q);

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
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
        reactionProbabilities[2] = 1.0;
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
            probNorm[i] =
                reactionProbabilities[i]/totalReactionProbability;

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
                        // Atom ionisation is to occur
                        atomIonisationReaction = true;
                        break;
                    }
                    if (i == 2)
                    {
                        // Dissociation reaction is to occur
                        dissocReaction = true;
                        break;
                    }
                }
            }
        }
    }

    if (dissocReaction)
    {
        ++nDissociationReactionsPerTimeStep_;
        ++nTotReactionsDiss_;

        if (allowSplitting_)
        {
            relax_ = false;
            dissociateP(heatOfReactionDiss_, productIdsDiss_, p, q);
        }
    }
    if (ionisationReaction)
    {
        ++nIonisationReactionsPerTimeStep_;
        ++nTotReactionsIon_;

        if (allowSplitting_)
        {
            relax_ = false;
            ioniseP(heatOfReactionIon_, productIdsIon_, p, q);
        }
    }
    if (atomIonisationReaction)
    {
        ++nIonisationReactions2PerTimeStep_;
        ++nTotReactionsIon2_;

        if (allowSplitting_)
        {
            relax_ = false;
            ioniseQ(heatOfReactionIon2_, productIdsIon2_, p, q);
        }
    }
}


bool Foam::moleculeAtomDissociationIonisation::outputResults
(
    const label counterIndex
)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIdsDiss_[0]];
    const word& productB = cloud_.typeIdList()[productIdsDiss_[1]];

    const word& productC = cloud_.typeIdList()[productIdsIon_[0]];
    const word& productD = cloud_.typeIdList()[productIdsIon_[1]];

    const word& productE = cloud_.typeIdList()[productIdsIon2_[0]];
    const word& productF = cloud_.typeIdList()[productIdsIon2_[1]];

    if (write)
    {
        const List<DynamicList<uniGasParcel*>>& cellOccupancy
            = cloud_.cellOccupancy();

        List<label> mols(2, 0);
        volume_ = 0.0;

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
        label nTotReactionsDiss = nTotReactionsDiss_;
        label nTotReactionsIon = nTotReactionsIon_;
        label nTotReactionsIon2 = nTotReactionsIon2_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateDiss = 0.0;
            scalar reactionRateIon = 0.0;
            scalar reactionRateIon2 = 0.0;

            reactionRateDiss =
            (
                nTotReactionsDiss
               *cloud_.nParticle()
            )
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Dissociation type II reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + " << reactantB
                << ", reaction rate = " << reactionRateDiss
                << endl;

            reactionRateIon =
            (nTotReactionsIon*cloud_.nParticle())
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productC << " + " << productD << " + " << reactantB
                << ", reaction rate = " << reactionRateIon
                << endl;

            reactionRateIon2 =
                (nTotReactionsIon2*cloud_.nParticle())
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << reactantA << " + " << productE << " + " << productF
                << ", reaction rate = " << reactionRateIon2
                << endl;
        }
    }
    else
    {
        label nTotReactionsDiss = nTotReactionsDiss_;
        label nTotReactionsIon = nTotReactionsIon_;
        label nTotReactionsIon2 = nTotReactionsIon2_;

        label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep =
            nIonisationReactions2PerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());

            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }

        if (nTotReactionsDiss > VSMALL)
        {
            Info<< "Dissociation type II reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIon > VSMALL)
        {
            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIon2 > VSMALL)
        {
            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << reactantA << " + " << productA << " + "
                << productB
                << " is active, nReactions this time step = "
                << nIonisationReactions2PerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
