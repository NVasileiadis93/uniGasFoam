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

#include "atomAtomIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(atomAtomIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    atomAtomIonisationDissimilarSpecies,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atomAtomIonisationDissimilarSpecies::atomAtomIonisationDissimilarSpecies
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    productIdsIon_(),
    chargedAtom_(false),
    heatOfReactionIon_(),
    heatOfReactionIon2_(),
    nTotIonisationReactions_(0),
    nTotIonisationReactions2_(0),
    nIonisationReactionsPerTimeStep_(0),
    nIonisationReactions2PerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atomAtomIonisationDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::atomAtomIonisationDissimilarSpecies::setProperties()
{
    if (reactants_[0] == reactants_[1])
    {
        FatalErrorInFunction
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    chargedAtom_ = propsDict_.get<bool>("chargedAtom");

    // check that reactant one is an 'atom'

    if (rDof1_ > VSMALL)
    {
        FatalErrorInFunction
            << "First reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[0]
            << exit(FatalError);
    }

    if (vDof1_ > VSMALL)
    {
         FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    // check that reactant two is an 'atom'

    if (rDof2_ > VSMALL)
    {
        FatalErrorInFunction
            << "Second reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[1]
            << exit(FatalError);
    }

    if (vDof2_ > VSMALL)
    {
         FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[1]
            << exit(FatalError);
    }

    //check that reactant one is not an 'ion' or an 'electron'

    if (charge1_ == -1)
    {
        FatalErrorInFunction
            << "First reactant cannot be an ion or an electron"
            << reactants_[0]
            << exit(FatalError);
    }

    //check that reactant two is not an 'ion' or an 'electron'

    if (charge2_ == -1)
    {
        FatalErrorInFunction
            << "Second reactant cannot be an ion or an electron"
            << reactants_[1]
            << exit(FatalError);
    }

    if (chargedAtom_)
    {
        // reading in ionisation products

        wordPair productsOfIonisedAtom
        (
            propsDict_.lookup("productsOfIonisedAtom")
        );

        productIdsIon_ = labelPair(-1, -1);

        forAll(productsOfIonisedAtom, r)
        {
            forAll(productIdsIon_, r)
            {
                productIdsIon_[r] =
                    cloud_.typeIdList().find(productsOfIonisedAtom[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon_[r] == -1)
                {
                    FatalErrorInFunction
                        << "Cannot find type id: "
                        << productsOfIonisedAtom[r]
                        << exit(FatalError);
                }
            }

            // check that product three is a 'atom'

            label rDof3 = cloud_.constProps(productIdsIon_[0]).rotationalDoF();

            label vDof3 = cloud_.constProps(productIdsIon_[0]).vibrationalDoF();

            if (rDof3 > 1)
            {
                FatalErrorInFunction
                    << "First product must be an atom (not a molecule): "
                    << productsOfIonisedAtom[0]
                    << exit(FatalError);
            }

            if (vDof3 > VSMALL)
            {
                FatalErrorInFunction
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:" << productIdsIon_[1]
                    << exit(FatalError);
            }

            // check that product two is an 'electron'

            const label charge = cloud_.constProps(productIdsIon_[1]).charge();

            if (charge != -1)
            {
                FatalErrorInFunction
                    << "Second product must be an electron: "
                    << productsOfIonisedAtom[1]
                    << exit(FatalError);
            }
        }
    }
    else
    {
        // reading in ionisation products

        wordPair productsOfIonisedAtom
        (
            propsDict_.lookup("productsOfIonisedAtomP")
        );

        productIdsIon_ = labelPair(-1, -1);

        forAll(productsOfIonisedAtom, r)
        {
            forAll(productIdsIon_, r)
            {
                productIdsIon_[r] =
                    cloud_.typeIdList().find(productsOfIonisedAtom[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon_[r] == -1)
                {
                    FatalErrorInFunction
                        << "Cannot find type id: "
                        << productsOfIonisedAtom[r]
                        << exit(FatalError);
                }
            }

            // check that product one is a 'atom'

            label rDof4 = cloud_.constProps(productIdsIon_[0]).rotationalDoF();

            label vDof4 =
                cloud_.constProps(productIdsIon_[0]).vibrationalDoF();

            if (rDof4 > 1)
            {
                FatalErrorInFunction
                    << "First product must be an atom (not an atom): "
                    << productsOfIonisedAtom[0]
                    << exit(FatalError);
            }

            if (vDof4 > VSMALL)
            {
                FatalErrorInFunction
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:"
                    << productsOfIonisedAtom[0]
                    << exit(FatalError);
            }

            // check that product two is an 'electron'

            label charge = cloud_.constProps(productIdsIon_[1]).charge();

            if (charge != -1)
            {
                FatalErrorInFunction
                    << "Second product must be an electron: "
                    << productsOfIonisedAtom[1]
                    << exit(FatalError);
            }
        }

        wordPair productsOfIonisedAtom2
        (
            propsDict_.lookup("productsOfIonisedAtomQ")
        );

        productIdsIon2_.setSize(productsOfIonisedAtom2.size());

        forAll(productsOfIonisedAtom2, r)
        {
            forAll(productIdsIon2_, r)
            {
                productIdsIon2_[r] =
                    cloud_.typeIdList().find(productsOfIonisedAtom2[r]);

                // check that reactants belong to the typeIdList
                if (productIdsIon2_[r] == -1)
                {
                    FatalErrorInFunction
                        << "Cannot find type id: "
                        << productsOfIonisedAtom2[r]
                        << exit(FatalError);
                }
            }

            // check that product one is a 'atom', not an 'molecule'

            label rDof5 =
                cloud_.constProps(productIdsIon2_[0]).rotationalDoF();

            label vDof5 =
                cloud_.constProps(productIdsIon2_[0]).vibrationalDoF();

            if (rDof5 > 1)
            {
                FatalErrorInFunction
                    << "First product must be an atom (not an atom): "
                    << productsOfIonisedAtom2[0]
                    << exit(FatalError);
            }

            if (vDof5 > VSMALL)
            {
                FatalErrorInFunction
                    << "Reactions are currently only implemented for "
                    << "monatomic and diatomic species"
                    << " This is a polyatomic:"
                    << productsOfIonisedAtom2[0]
                    << exit(FatalError);
            }

            // check that product two is an 'electron'

            const label charge = cloud_.constProps(productIdsIon2_[1]).charge();

            if (charge != -1)
            {
                FatalErrorInFunction
                    << "Second product must be an electron: "
                    << productsOfIonisedAtom2[1]
                    << exit(FatalError);
            }
        }
    }

    if (chargedAtom_)
    {
        heatOfReactionIon_ =
            propsDict_.get<scalar>("heatOfReactionIonisation");
    }
    else
    {
        heatOfReactionIon_ =
            propsDict_.get<scalar>("heatOfReactionIonisationP");
        heatOfReactionIon2_ =
            propsDict_.get<scalar>("heatOfReactionIonisationQ");
    }
}


bool Foam::atomAtomIonisationDissimilarSpecies::tryReactMolecules
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


void Foam::atomAtomIonisationDissimilarSpecies::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    if (chargedAtom_)
    {
        if (p.typeId() == reactantIds_[0] && p.typeId() == reactantIds_[1])
        {
            reactCharged
            (
                p,
                q,
                heatOfReactionIon_,
                productIdsIon_,
                nTotIonisationReactions_,
                nIonisationReactionsPerTimeStep_
            );
        }
        else if (p.typeId() == reactantIds_[1] && q.typeId() == reactantIds_[0])
        {
            reactCharged
            (
                q,
                p,
                heatOfReactionIon_,
                productIdsIon_,
                nTotIonisationReactions_,
                nIonisationReactionsPerTimeStep_
            );
        }
    }
    else
    {
        if (p.typeId() == reactantIds_[0] && q.typeId() == reactantIds_[1])
        {
            reactNotCharged
            (
                p,
                q,
                heatOfReactionIon_,
                productIdsIon_,
                heatOfReactionIon2_,
                productIdsIon2_,
                nTotIonisationReactions_,
                nIonisationReactionsPerTimeStep_,
                nTotIonisationReactions2_,
                nIonisationReactions2PerTimeStep_
            );
        }
        else if (p.typeId() == reactantIds_[1] && q.typeId() == reactantIds_[0])
        {
            reactNotCharged
            (
                q,
                p,
                heatOfReactionIon2_,
                productIdsIon2_,
                heatOfReactionIon_,
                productIdsIon_,
                nTotIonisationReactions2_,
                nIonisationReactions2PerTimeStep_,
                nTotIonisationReactions_,
                nIonisationReactionsPerTimeStep_
            );
        }
    }
}


void Foam::atomAtomIonisationDissimilarSpecies::reactNotCharged
(
    uniGasParcel& p,
    uniGasParcel& q,
    const scalar heatOfReaction,
    const labelPair& productIds,
    const scalar heatOfReaction2,
    const labelPair& productIds2,
    label& nReactions,
    label& nReactionsPerTimeStep,
    label& nReactions2,
    label& nReactions2PerTimeStep
)
{
    relax_ = true;

    scalar totalProbability = 0.0;
    FixedList<scalar, 2> reactionProbabilities(0.0);

    scalar EcPPIon = 0.0;
    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
       *physicoChemical::k.value();

    // Calculate if an ionisation of P is possible
    EcPPIon = translationalEnergy(p, q) + EEle(p);

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        totalProbability += 1.0;
        reactionProbabilities[0] = 1.0;
    }

    ionisationEnergy =
        cloud_.constProps(q.typeId()).ionisationTemperature()
       *physicoChemical::k.value();

    // Calculate if an ionisation of Q is possible
    EcPPIon = translationalEnergy(p, q) + EEle(q);

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        totalProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }

    // Decide if a reaction is to occur

    if (totalProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        FixedList<scalar, 2> probNorm;
        probNorm[0] = reactionProbabilities[0]/totalProbability;
        probNorm[1] = reactionProbabilities[1]/totalProbability;

        scalar sumProbability = 0.0;

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
                        ++nReactions;
                        ++nReactionsPerTimeStep;

                        if (allowSplitting_)
                        {
                            relax_ = false;
                            ioniseP(heatOfReaction, productIds, p, q);
                        }

                        break;
                    }
                    if (i == 1)
                    {
                        ++nReactions2;
                        ++nReactions2PerTimeStep;

                        if (allowSplitting_)
                        {
                            relax_ = false;
                            ioniseQ(heatOfReaction, productIds, p, q);
                        }

                        break;
                    }
                }
            }
        }
    }
}


void Foam::atomAtomIonisationDissimilarSpecies::reactCharged
(
    uniGasParcel& p,
    uniGasParcel& q,
    const scalar heatOfReaction,
    const labelPair& productIds,
    label& nReactions,
    label& nReactionsPerTimeStep
)
{
    relax_ = true;

    // Calculate if an ionisation of P is possible
    scalar EcPPIon = translationalEnergy(p, q) + EEle(p);

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
       *physicoChemical::k.value();

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        ++nReactions;
        ++nReactionsPerTimeStep;

        if (allowSplitting_)
        {
            relax_ = false;
            ioniseP(heatOfReaction, productIds, p, q);
        }
    }
}


bool Foam::atomAtomIonisationDissimilarSpecies::outputResults
(
    const label counterIndex
)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
    
    const word& productA = cloud_.typeIdList()[productIdsIon_[0]];
    const word& productB = cloud_.typeIdList()[productIdsIon_[1]];
    
    word productC;
    word productD;
    
    if (!chargedAtom_)
    {
        productC = cloud_.typeIdList()[productIdsIon2_[0]];
        productD = cloud_.typeIdList()[productIdsIon2_[1]];
    }

//     const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
//     const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
// 
//     const word& productA = cloud_.typeIdList()[productIdsIon_[0]];
//     const word& productB = cloud_.typeIdList()[productIdsIon_[1]];
// 
//     word productC;
//     word productD;

    if (!chargedAtom_)
    {
        productC = cloud_.typeIdList()[productIdsIon2_[0]];
        productD = cloud_.typeIdList()[productIdsIon2_[1]];
    }

    if (write)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        List<label> mols(2, 0);
        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            for (uniGasParcel* p : cellOccupancy[c])
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
        label nTotReactionsIonisation = nTotIonisationReactions_;
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0) && (numberDensities_[1] > 0))
        {
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateIonisation2 = 0.0;

            reactionRateIonisation =
                (nTotReactionsIonisation*cloud_.nParticle())
               /(
                    counterIndex
                   *deltaT
                   *numberDensities_[0]
                   *numberDensities_[1]
                   *volume
                );

            Info<< "Ionisation reaction "
                << reactantA << " + " << reactantB
                << " --> "
                << productA << " + " << productB << " + " << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            if (!chargedAtom_)
            {
                reactionRateIonisation2 =
                    (nTotReactionsIonisation2*cloud_.nParticle())
                   /(
                        counterIndex
                       *deltaT
                       *numberDensities_[0]
                       *numberDensities_[1]
                       *volume
                    );

                Info<< "Ionisation reaction "
                    << reactantA << " + " << reactantB
                    << " --> "
                    << reactantA << " + " << productC
                    << " + " << productD
                    << ", reaction rate = " << reactionRateIonisation2
                    << endl;
            }
        }
    }
    else
    {
        label nTotReactionsIonisation = nTotIonisationReactions_;
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep =
            nIonisationReactions2PerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }

        if (nTotReactionsIonisation > VSMALL)
        {
            Info<< "Ionisation reaction "
                << reactantA << " + " << reactantB
                << " --> "
                << productA << " + " << productB
                << " + " << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }

        if (nTotReactionsIonisation2 > VSMALL)
        {
            Info<< "Ionisation reaction "
                << reactantA << " + " << reactantB
                << " --> "
                << reactantA << " + " << productC
                << " + " << productD
                << " is active, nReactions this time step = "
                << nIonisationReactions2PerTimeStep << endl;
        }
    }

    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
