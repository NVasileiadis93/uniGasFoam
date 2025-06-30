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

#include "chargeExchange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(chargeExchange, 0);
addToRunTimeSelectionTable(uniGasReaction, chargeExchange, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chargeExchange::chargeExchange
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductIds_(),
    dissociationPossible_(propsDict_.get<bool>("dissociationPossible")),
    heatOfReactionDissoc_(0.0),
    heatOfReactionIon_(propsDict_.get<scalar>("heatOfReactionIonisation")),
    heatOfReactionExch_(propsDict_.get<scalar>("heatOfReactionExchange")),
    activationEnergy_(0.0),
    aCoeff_(propsDict_.get<scalar>("aCoeff")),
    bCoeff_(propsDict_.get<scalar>("bCoeff")),
    nChargeExchangeReactions_(0),
    nIonisationReactions_(0),
    nDissociationReactions_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chargeExchange::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::chargeExchange::setProperties()
{
    // Check that reactants are 'atoms' or 'molecules', not 'electrons'

    if (charge1_ == -1)
    {
        FatalErrorInFunction
            << "Reactant must not be an electron: " << reactants_[0]
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

    if (charge2_ == -1)
    {
        FatalErrorInFunction
            << "Reactant must not be an electron: " << reactants_[1]
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

    // Reading in charge exchange products

    wordPair chargeExchangeProductMolecules
    (
        propsDict_.lookup("chargeExchangeProducts")
    );

    chargeExchangeProductIds_ = labelPair(-1, -1);

    forAll(chargeExchangeProductIds_, i)
    {
        chargeExchangeProductIds_[i] =
            cloud_.typeIdList().find(chargeExchangeProductMolecules[i]);

        // Check that products belong to the typeIdList
        if (chargeExchangeProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: "
                << chargeExchangeProductMolecules[i] << nl
                << exit(FatalError);
        }

        // Check that products are a 'molecule' or an 'atom'

        label charge = cloud_.constProps(chargeExchangeProductIds_[i]).charge();

        if (charge == -1)
        {
            FatalErrorInFunction
                << "Products cannot be an electron: "
                << chargeExchangeProductMolecules
                << exit(FatalError);
        }
    }

    // Reading in dissociation products

    if (dissociationPossible_)
    {
        const wordPair dissociationProducts
        (
            propsDict_.lookup("dissociationProducts")
        );

        dissociationProductIds_ = labelPair(-1, -1);

        forAll(dissociationProductIds_, i)
        {
            dissociationProductIds_[i] =
                cloud_.typeIdList().find(dissociationProducts[i]);

            // check that products belong to the typeIdList
            if (dissociationProductIds_[i] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << dissociationProducts[i] << nl
                    << exit(FatalError);
            }

            // check that products are 'atoms'

            label iD = dissociationProductIds_[i];

            label rDof = cloud_.constProps(iD).rotationalDoF();

            if (rDof > VSMALL)
            {
                FatalErrorInFunction
                    << "Dissociation products must be atoms: "
                    << dissociationProducts
                    << exit(FatalError);
            }
        }
    }

    // Reading in ionisation products

    wordPair ionisationProductMolecules
    (
        propsDict_.lookup("ionisationProducts")
    );

    ionisationProductIds_ = labelPair(-1, -1);

    forAll(ionisationProductIds_, i)
    {
        ionisationProductIds_[i] =
            cloud_.typeIdList().find(ionisationProductMolecules[i]);

        // Check that products belong to the typeIdList
        if (ionisationProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: "
                << ionisationProductMolecules[i] << nl
                << exit(FatalError);
        }
    }

    // Check that second product is an 'electron'

    const label charge = cloud_.constProps(ionisationProductIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second ionisation product must be an electron: "
            << ionisationProductMolecules
            << exit(FatalError);
    }

    activationEnergy_ *= physicoChemical::k.value();

    if (dissociationPossible_)
    {
        heatOfReactionDissoc_ =
            propsDict_.get<scalar>("heatOfReactionDissociation");
    }
}


bool Foam::chargeExchange::tryReactMolecules
(
    const label typeIdP,
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


void Foam::chargeExchange::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    if (p.typeId() == reactantIds_[0] && p.typeId() == reactantIds_[1])
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

void Foam::chargeExchange::react
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
    bool chargeExchange = false;

    // 3 reactions possible:
    // - 1. Dissociation
    // - 2. Ionisation
    // - 3. Charge exchange

    scalar EcP = 0.0;
    scalar TColl = 0.0;

    // Check for dissociation

    if (dissociationPossible_)
    {
        if (cloud_.constProps(p.typeId()).rotationalDoF() > 0)
        {
            label idP = charDissLevel(p);
            label imaxP = 0;

            EcP = translationalEnergy(p, q) + EVib(p);

            imaxP = EcP/(physicoChemical::k.value()*thetaV(p));

            if (imaxP - idP > 0)
            {
                // Dissociation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[0] = 1.0;
            }
        }
    }

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation is possible
        
    EcP = translationalEnergy(p, q) + EEle(p);

    if ((EcP - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }
    
    // calculate the charge exchange probability

    scalar heatOfReactionExchJoules =
        heatOfReactionExch_*physicoChemical::k.value();

    scalar aDash =
        aCoeff_*(pow(2.5 - omega(p, q), bCoeff_)
        *exp(lgamma(2.5 - omega(p, q)))
        /exp(lgamma(2.5 - omega(p, q) + bCoeff_)));

    TColl =
        (translationalEnergy(p, q)/(physicoChemical::k.value()))
        /(2.5 - omega(p, q));

    scalar activEn =
        activationEnergy_
        + (aDash*pow((TColl/273.0), bCoeff_)
        * mag(heatOfReactionExchJoules));

    if (heatOfReactionExchJoules < 0.0)
    {
        activEn -= heatOfReactionExchJoules;
    }

    if (EcP > activEn)
    {
        label keyElectronicLevel = -1;

        for (label i = 0; i < jMax(p) + 1; ++i)
        {
            scalar electronicEnergy = EEList(p)[i];

            if (electronicEnergy > activEn)
            {
                break;
            }

            ++keyElectronicLevel;
        }

        label trialELevel =
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcP,
                jMax(p) + 1,
                omega(p, q),
                EEList(p),
                gList(p)
            );

        if (trialELevel == keyElectronicLevel)
        {
            scalar prob = 0.0;

            label nPossStates = 0;

            if ((jMax(p) + 1) == 1)
            {
                nPossStates = gList(p)[0];
            }
            else
            {
                forAll(EEList(p), n)
                {
                    if (EcP > EEList(p)[n])
                    {
                        nPossStates += gList(p)[n];
                    }
                }
            }

            label nState =
                ceil(cloud_.rndGen().sample01<scalar>()*nPossStates);
            label nAvailableStates = 0;
            label nLevel = -1;

            forAll(EEList(p), n)
            {
                nAvailableStates += gList(p)[n];

                if (nState <= nAvailableStates && nLevel < 0)
                {
                    nLevel = n;
                }
            }

            // Calculate the probability of it occurring
            scalar summation = 0.0;

            for (label i = 0; i <= nLevel; ++i)
            {
                summation +=
                    gList(p)[i]
                    *pow(EcP - EEList(p)[i], 1.5 - omega(p, q));
            }

            prob =
                gList(p)[trialELevel]
                *pow
                (
                    EcP - EEList(p)[trialELevel],
                    1.5 - omega(p, q)
                )
                /summation;

            if (prob > cloud_.rndGen().sample01<scalar>())
            {
                // Charge exchange can occur
                totalReactionProbability += prob;
                reactionProbabilities[2] = prob;
            }
        }
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm(reactionProbabilities/totalReactionProbability);
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
                        // Dissociation is to occur
                        dissocReaction = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Ionisation reaction is to occur
                        ionisationReaction = true;
                        break;
                    }
                    if (i == 2)
                    {
                        // Charge exchange is to occur
                        chargeExchange = true;
                        break;
                    }
                }
            }
        }
    }

    if (dissocReaction)
    {
        ++nDissociationReactions_;
        ++nDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDissoc_,
                dissociationProductIds_,
                p,
                q
            );
        }
    }

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
                ionisationProductIds_,
                p,
                q
            );
        }
    }

    if (chargeExchange)
    {
        ++nChargeExchangeReactions_;
        ++nChargeExchangeReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;
            
            chargeExchangePQ
            (
                heatOfReactionExchJoules,
                chargeExchangeProductIds_,
                p,
                q
            );
        }
    }
}


bool Foam::chargeExchange::outputResults(const label counterIndex)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    word reactantA = cloud_.typeIdList()[reactantIds_[0]];
    word reactantB = cloud_.typeIdList()[reactantIds_[1]];

    word productA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
    word productB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];

    word productC;
    word productD;

    if (dissociationPossible_)
    {
        productC = cloud_.typeIdList()[dissociationProductIds_[0]];
        productD = cloud_.typeIdList()[dissociationProductIds_[1]];
    }

    word productE = cloud_.typeIdList()[ionisationProductIds_[0]];
    word productF = cloud_.typeIdList()[ionisationProductIds_[1]];
    
    label nTotChargeExchangeReactions = nChargeExchangeReactions_;
    label nTotDissociationReactions = nDissociationReactions_;
    label nTotIonisationReactions = nIonisationReactions_;

    label nChargeExchangeReactionsPerTimeStep =
        nChargeExchangeReactionsPerTimeStep_;
    label nDissociationReactionsPerTimeStep =
        nDissociationReactionsPerTimeStep_;
    label nIonisationReactionsPerTimeStep =
        nIonisationReactionsPerTimeStep_;

    if (Pstream::parRun())
    {
        reduce(nTotChargeExchangeReactions, sumOp<label>());
        reduce(nTotDissociationReactions, sumOp<label>());
        reduce(nTotIonisationReactions, sumOp<label>());
        reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
        reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
        reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
    }

    if (write)
    {
        // measure density

        const auto& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols(2, 0);
        scalar volume = volume_;

        forAll(cellOccupancy, c)
        {
            const auto& parcelsInCell = cellOccupancy[c];

            for (uniGasParcel* p : parcelsInCell)
            {
                label id = reactantIds_.find(p->typeId());

                if (id != -1)
                {
                    mols[id]++;
                }
            }
        }

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(volume, sumOp<scalar>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;

        if (reactantIds_[0] == reactantIds_[1])
        {
            numberDensities_[1] = (mols[0]*cloud().nParticle())/volume;
        }
        else
        {
            numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;
        }

        const scalar deltaT = mesh_.time().deltaTValue();


        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateChargeExchange = 0.0;
            scalar reactionRateDissociation = 0.0;
            scalar reactionRateIonisation = 0.0;

            reactionRateChargeExchange =
                (nTotChargeExchangeReactions*cloud_.nParticle())
               /(
                    counterIndex
                   *deltaT
                   *numberDensities_[0]
                   *numberDensities_[1]
                   *volume
                );

            Info<< "Charge exchange reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB
                << ", reaction rate = " << reactionRateChargeExchange
                << endl;

            reactionRateIonisation =
                (nTotIonisationReactions*cloud_.nParticle())
               /(
                    counterIndex
                   *deltaT
                   *numberDensities_[0]
                   *numberDensities_[1]
                   *volume
                );

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productE << " + " << productF
                 << " + " << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;

            if (dissociationPossible_)
            {
                reactionRateDissociation =
                    (nTotDissociationReactions*cloud_.nParticle())
                   /(
                        counterIndex
                       *deltaT
                       *numberDensities_[0]
                       *numberDensities_[1]
                       *volume
                    );

                Info<< "Dissociation reaction "
                    <<  reactantA << " + " << reactantB
                    <<  " --> "
                    << productC << " + " << productD
                    << " + " << reactantB
                    << ", reaction rate = " << reactionRateDissociation
                    << endl;
            }
        }
    }
    else
    {
        if (nTotChargeExchangeReactions > VSMALL)
        {
                    Info<< "Charge exchange reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB
                << " is active, nReactions this time step = "
                << nChargeExchangeReactionsPerTimeStep << endl;
        }

        if (nTotDissociationReactions > VSMALL)
        {
            Info<< "Dissociation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productC << " + " << productD
                << " + " << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactions > VSMALL)
        {

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productE << " + " << productF
                << " + " << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nChargeExchangeReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
