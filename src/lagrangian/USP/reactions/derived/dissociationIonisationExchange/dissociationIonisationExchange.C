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

#include "dissociationIonisationExchange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(dissociationIonisationExchange, 0);

addToRunTimeSelectionTable
(
    uspReaction,
    dissociationIonisationExchange,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dissociationIonisationExchange::dissociationIonisationExchange
(
    const Time& t,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspReaction(t, cloud, dict),
    exchangeProductIds_(),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductsIdsP_(),
    ionisationProductsIdsQ_(),
    chargedAtom_(propsDict_.get<bool>("chargedAtom")),
    chargedMolecule_(propsDict_.get<bool>("chargedMolecule")),
    chargeExchange_(propsDict_.get<bool>("chargeExchange")),
    heatOfReactionDiss_(),
    heatOfReactionExch_(propsDict_.get<scalar>("heatOfReactionExch")),
    heatOfReactionChargeExchange_(),
    heatOfReactionIonP_(),
    heatOfReactionIonQ_(),
    aCoeff_(propsDict_.get<scalar>("aCoeff")),
    bCoeff_(propsDict_.get<scalar>("bCoeff")),
    aCoeffCharge_(),
    bCoeffCharge_(),
    nTotExchangeReactions_(0),
    nTotChargeExchangeReactions_(),
    nTotDissociationReactions_(0),
    nTotIonisationReactionsP_(0),
    nTotIonisationReactionsQ_(0),
    nExchangeReactionsPerTimeStep_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPPerTimeStep_(0),
    nIonisationReactionsQPerTimeStep_(0)
{
    if (!chargedMolecule_)
    {
        heatOfReactionDiss_ = propsDict_.get<scalar>("heatOfReactionDiss");
        heatOfReactionIonP_ = propsDict_.get<scalar>("heatOfReactionIonP");
    }

    if (!chargedAtom_)
    {
        heatOfReactionIonQ_ = propsDict_.get<scalar>("heatOfReactionIonQ");
    }

    if (chargeExchange_)
    {
        heatOfReactionChargeExchange_ =
            propsDict_.get<scalar>("heatOfReactionChargeExch");
        aCoeffCharge_ = propsDict_.get<scalar>("aCoeffChargeExchange");
        bCoeffCharge_ = propsDict_.get<scalar>("bCoeffChargeExchange");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dissociationIonisationExchange::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::dissociationIonisationExchange::setProperties()
{
    // check that the first reactant is a 'molecule'

    if (rDof1_ < VSMALL)
    {
        FatalErrorInFunction
            << "Reactant 1 must be a molecule (not an atom): "
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

    // check that the second reactant is an 'atom'

    if (rDof2_ > VSMALL)
    {
        FatalErrorInFunction
            << "Reactant 2 must be an atom (not a molecule): "
            << reactants_[1]
            << exit(FatalError);
    }

    //reading in products

    const wordPair exchangeProductMolecules
    (
        propsDict_.lookup("productsOfExchangeReaction")
    );

    if (exchangeProductMolecules[0] == exchangeProductMolecules[1])
    {
        FatalErrorInFunction
            << "Exchange reaction product molecules cannot be same species."
            << exit(FatalError);
    }

    exchangeProductIds_ = labelPair(-1, -1);

    forAll(exchangeProductIds_, r)
    {
        exchangeProductIds_[r] =
            cloud_.typeIdList().find(exchangeProductMolecules[r]);

        // check that reactants belong to the typeIdList
        if (exchangeProductIds_[r] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << exchangeProductMolecules[r]
                << exit(FatalError);
        }
    }

    // check that first exchange product is a 'molecule'

    const label rDof3 =
        cloud_.constProps(exchangeProductIds_[0]).rotationalDoF();

    const label vDof3 =
        cloud_.constProps(exchangeProductIds_[0]).vibrationalDoF();

    if (rDof3 < 0)
    {
        FatalErrorInFunction
            << "First product of the exchange reaction must"
            << " be a molecule (not an atom): "
            << exchangeProductMolecules[0]
            << exit(FatalError);
    }

    if (vDof3 > 1)
    {
         FatalErrorInFunction
            << "Reactions are currently only implemented for "
            << "monatomic and diatomic species."
            << " This is a polyatomic:" << exchangeProductMolecules[0]
            << exit(FatalError);
    }

    // check that second exchange product is an 'atom'

    label id = exchangeProductIds_[1];

    const label rDof4 = cloud_.constProps(id).rotationalDoF();

    if (rDof4 > 0)
    {
        FatalErrorInFunction
            << "Second product of the exchange reaction "
            << "must be an atom (not a molecule): "
            << exchangeProductMolecules[1]
            << exit(FatalError);
    }

    if (chargeExchange_)
    {
        // reading in products

        const wordPair chargeExchangeProductMolecules
        (
            propsDict_.lookup("productsOfChargeExchangeReaction")
        );

        if
        (
            chargeExchangeProductMolecules[0]
         == chargeExchangeProductMolecules[1]
        )
        {
            FatalErrorInFunction
                << "Charge exchange reaction product molecules "
                << "cannot be same species."
                << exit(FatalError);
        }

        chargeExchangeProductIds_ = labelPair(-1, -1);

        forAll(chargeExchangeProductIds_, r)
        {
            chargeExchangeProductIds_[r] =
                cloud_.typeIdList().find(chargeExchangeProductMolecules[r]);

            // check that reactants belong to the typeIdList
            if (chargeExchangeProductIds_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << chargeExchangeProductMolecules[r]
                    << exit(FatalError);
            }
        }

        // check that first exchange product is a 'molecule'

        label id2 = chargeExchangeProductIds_[0];

        const label rDof5 = cloud_.constProps(id2).rotationalDoF();

        if (rDof5 < VSMALL)
        {
            FatalErrorInFunction
                << "First product of the exchange reaction must "
                << "be a molecule (not an atom): "
                << chargeExchangeProductMolecules[0]
                << exit(FatalError);
        }

        // check that second exchange product is an 'atom'

        label id3 = chargeExchangeProductIds_[1];

        const label rDof6 = cloud_.constProps(id3).rotationalDoF();

        if (rDof6 > VSMALL)
        {
            FatalErrorInFunction
                << "Second product of the exchange reaction "
                << "must be an atom (not a molecule): "
                << chargeExchangeProductMolecules[1]
                << exit(FatalError);
        }
    }

    if (!chargedMolecule_)
    {
        const wordPair dissociationProductMolecules
        (
            propsDict_.lookup("productsOfDissociatedMolecule")
        );

        dissociationProductIds_ = labelPair(-1, -1);

        forAll(dissociationProductIds_, r)
        {
            dissociationProductIds_[r] =
                cloud_.typeIdList().find(dissociationProductMolecules[r]);

            // check that reactants belong to the typeIdList
            if (dissociationProductIds_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << dissociationProductMolecules[r]
                    << exit(FatalError);
            }
        }

        label id4 = dissociationProductIds_[0];

        const label rDof7 = cloud_.constProps(id4).rotationalDoF();

        if (rDof7 > VSMALL)
        {
            FatalErrorInFunction
                << "First product of the dissociation reaction must "
                << "be an atom (not a molecule): "
                << dissociationProductMolecules[0]
                << exit(FatalError);
        }

        // check that second exchange product is an 'atom'

        label id5 = dissociationProductIds_[1];

        const label rDof8 = cloud_.constProps(id5).rotationalDoF();

        if (rDof8 > VSMALL)
        {
            FatalErrorInFunction
                << "Second product of the exchange reaction must "
                << "be an atom (not a molecule): "
                << dissociationProductMolecules[1]
                << exit(FatalError);
        }

        // read in ionisation products

        const wordPair ionisationProductMoleculesP
        (
            propsDict_.lookup("productsOfIonisedMolecule")
        );

        ionisationProductsIdsP_ = labelPair(-1, -1);

        forAll(ionisationProductsIdsP_, r)
        {
            ionisationProductsIdsP_[r] =
                cloud_.typeIdList().find(ionisationProductMoleculesP[r]);

            // check that reactants belong to the typeIdList
            if (ionisationProductsIdsP_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << ionisationProductMoleculesP[r]
                    << exit(FatalError);
            }
        }

        // check that first ionisation product is a 'molecule'

        label id6 = ionisationProductsIdsP_[0];

        const label rDof9 = cloud_.constProps(id6).rotationalDoF();

        if (rDof9 < VSMALL)
        {
            FatalErrorInFunction
                << "First product of the molecule ionisation reaction "
                << "must be a molecule (not an atom): "
                << ionisationProductMoleculesP[0]
                << exit(FatalError);
        }

        // check that second ionisation product is a 'electron'

        label id7 = ionisationProductsIdsP_[1];

        const label charge = cloud_.constProps(id7).charge();

        if (charge != -1)
        {
            FatalErrorInFunction
                << "Second product of the molecule ionisation reaction must "
                << "be an electron: " << ionisationProductMoleculesP[1]
                << exit(FatalError);
        }
    }

    if (!chargedAtom_)
    {
        const wordPair ionisationProductMoleculesQ
        (
            propsDict_.lookup("productsOfIonisedAtom")
        );

        ionisationProductsIdsQ_ = labelPair(-1, -1);

        forAll(ionisationProductsIdsQ_, r)
        {
            ionisationProductsIdsQ_[r] =
                cloud_.typeIdList().find(ionisationProductMoleculesQ[r]);

            // check that reactants belong to the typeIdList
            if (ionisationProductsIdsQ_[r] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find type id: "
                    << ionisationProductMoleculesQ[r]
                    << nl
                    << exit(FatalError);
            }
        }

        // check that first ionisation product is an 'atom'

        label id8 = ionisationProductsIdsQ_[0];

        const label rDof10 = cloud_.constProps(id8).rotationalDoF();

        if (rDof10 > VSMALL)
        {
            FatalErrorInFunction
                << "First product of the atom ionisation reaction must be "
                << "an atom (not a molecule): "
                << ionisationProductMoleculesQ[0]
                << exit(FatalError);
        }

        // check that second ionisation product is a 'electron'

        label id9 = ionisationProductsIdsQ_[1];

        const label charge = cloud_.constProps(id9).charge();

        if (charge != -1)
        {
            FatalErrorInFunction
                << "Second product of the atom ionisation reaction must be "
                << "an electron: " << ionisationProductMoleculesQ[1]
                << exit(FatalError);
        }
    }
}


bool Foam::dissociationIonisationExchange::tryReactMolecules
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


void Foam::dissociationIonisationExchange::reaction
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


void Foam::dissociationIonisationExchange::react
(
    uspParcel& p,
    uspParcel& q
)
{
    relax_ = true;

    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(5, 0.0);

    scalar EcPQ = 0.0;
    scalar TColl = 0.0;
    label idP = -1;
    label deltaDissoIP = 0;
    label imaxP = 0;
    label iaP = 0;
    bool dissocReaction = false;
    bool ionisationReactionP = false;
    bool ionisationReactionQ = false;
    bool exchangeReaction = false;
    bool chargeExchangeReaction = false;

    if (!chargedMolecule_)
    {
        // Calculate dissociation probability (0 or 1).
        EcPQ = translationalEnergy(p, q) + EVib(p);

        imaxP = EcPQ/(physicoChemical::k.value()*thetaV(p));

        idP = thetaD(p)/thetaV(p);

        deltaDissoIP = imaxP - idP;

        if (deltaDissoIP > 0)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        // Ionisation of the molecule

        scalar ionisationEnergy =
                        cloud_.constProps(p.typeId()).ionisationTemperature()
                        *physicoChemical::k.value();

        // calculate if an ionisation of the molecule is possible
        EcPQ = translationalEnergy(p, q) + EEle(p);

        if ((EcPQ - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
    }

    if (!chargedAtom_)
    {
        // Ionisation of the atom

        scalar ionisationEnergy =
                        cloud_.constProps(q.typeId()).ionisationTemperature()
                        *physicoChemical::k.value();

        // calculate if an ionisation of thee atom is possible
        EcPQ = translationalEnergy(p, q) + EEle(q);

        if ((EcPQ - ionisationEnergy) > VSMALL)
        {
            // Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }
    }

    EcPQ = translationalEnergy(p, q) + EVib(p);

    scalar P_exch = 0.0;

    // Next, calculate exchange probability.

    TColl =
        (translationalEnergy(p, q)/(physicoChemical::k.value()))
        /(2.5 - omega(p, q));

    scalar heatOfReactionExchJoules =
        heatOfReactionExch_*physicoChemical::k.value();

    scalar aDash =
        aCoeff_
        *(
            pow(2.5 - omega(p, q), bCoeff_)
            *exp(lgamma(2.5 - omega(p, q)))
            /exp(lgamma(2.5 - omega(p, q) + bCoeff_))
        );
        
    scalar activationEnergy = 0.0;

    if (heatOfReactionExchJoules < VSMALL)
    {
        //forward exchange
        activationEnergy = mag(heatOfReactionExchJoules)
                                *(1.0 + aDash*pow((TColl/273.0), bCoeff_));
    }
    else
    {
        //reverse exchange
        activationEnergy = (aDash*pow((TColl/273.0), bCoeff_)
                                * mag(heatOfReactionExchJoules));
    }

    if (EcPQ > activationEnergy)
    {
        scalar summation = 0.0;

        if
        (
            activationEnergy
            < cloud_.constProps(p.typeId()).thetaV()[0]
            *physicoChemical::k.value()
        )
        {
            summation = 1.0;
        }
        else
        {
            iaP = EcPQ/(physicoChemical::k.value()*thetaV(p));

            for (label i = 0; i <= iaP; ++i)
            {
                summation +=
                pow
                (
                    1.0 - (i*physicoChemical::k.value()*thetaV(p))/EcPQ,
                    1.5 - omega(p, q)
                );
            }
        }

        P_exch =
            pow
            (
                1.0 - (activationEnergy/EcPQ),
                1.5 - omega(p, q)
            )
            /summation;

        totalReactionProbability += P_exch;
        reactionProbabilities[3] = P_exch;
    }

    if (chargeExchange_)
    {
        // calculate charge exchange probability

        label maxElectronicLevelP =
            cloud_.constProps(p.typeId()).nElectronicLevels();

        scalar heatOfReactionChargeExchJoules =
            heatOfReactionChargeExchange_*physicoChemical::k.value();

        scalar EcP = translationalEnergy(p, q) + EEle(p);

        scalar aDash =
            aCoeffCharge_
            *(
                pow(2.5 - omega(p, q), bCoeffCharge_)
                *exp(lgamma(2.5 - omega(p, q)))
                /exp(lgamma(2.5 - omega(p, q) + bCoeffCharge_))
            );

        TColl =
            (translationalEnergy(p, q)/(physicoChemical::k.value()))
            /(2.5 - omega(p, q));

        scalar activationEnergy =
            aDash
            *pow((TColl/273.0), bCoeffCharge_)
            *mag(heatOfReactionChargeExchJoules);

        if (heatOfReactionChargeExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionChargeExchJoules;
        }

        if (EcP > activationEnergy)
        {
            label keyElectronicLevel = -1;

            for (label i = 0; i < maxElectronicLevelP; ++i)
            {
                scalar electronicEnergy = EEList(p)[i];

                if (electronicEnergy > activationEnergy)
                {
                    break;
                }

                ++keyElectronicLevel;
            }

            EcP = translationalEnergy(p, q) + EEle(p);

            label trialELevel =
                cloud_.postCollisionElectronicEnergyLevel
                (
                    EcP,
                    maxElectronicLevelP,
                    omega(p, q),
                    EEList(p),
                    gList(p)
                );

            if (trialELevel == keyElectronicLevel)
            {
                scalar prob = 0.0;

                label nPossStates = 0;

                if (maxElectronicLevelP == 1)
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

                // Calculate the probability of it occuring
                scalar summation = 0.0;

                for (label i = 0; i <= nLevel; i++)
                {
                    summation +=
                        gList(p)[i]
                        *pow(EcP - EEList(p)[i], 1.5 - omega(p, q));
                }

                prob =
                    (
                        gList(p)[trialELevel]
                        *pow
                        (
                            EcP - EEList(p)[trialELevel],
                            1.5 - omega(p, q))
                        )
                        /summation;

                if (prob > cloud_.rndGen().sample01<scalar>())
                {
                    // Charge exchange can occur
                    totalReactionProbability += prob;
                    reactionProbabilities[4] = prob;
                }
            }
        }

    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm
        (
            reactionProbabilities.size(),
            Zero
        );
        scalar sumProbability = 0.0;

        forAll(probNorm, i)
        {
            probNorm[i] =
                reactionProbabilities[i]/totalReactionProbability;

            // If current reaction can't occur, don't check for it
            if (probNorm[i] > VSMALL)
            {
                sumProbability += probNorm[i];

                if
                (
                    sumProbability
                    > cloud_.rndGen().sample01<scalar>()
                )
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
                        // Molecule ionisation reaction is to occur
                        ionisationReactionP = true;
                        break;
                    }
                    if (i == 2)
                    {
                        // Atom ionisation reaction is to occur
                        ionisationReactionQ = true;
                        break;
                    }
                    if (i == 3)
                    {
                        // Exchange reaction is to occur
                        exchangeReaction = true;
                        break;
                    }
                    if (i == 4)
                    {
                        // Exchange reaction is to occur
                        chargeExchangeReaction = true;
                        break;
                    }
                }
            }
        }
    }

    // Perform a dissociation reaction
    if (dissocReaction)
    {
        ++nTotDissociationReactions_;
        ++nDissociationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            dissociateP
            (
                heatOfReactionDiss_,
                dissociationProductIds_,
                p,
                q
            );
        }
    }

    if (ionisationReactionP)
    {
        ++nTotIonisationReactionsP_;
        ++nIonisationReactionsPPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonP_,
                ionisationProductsIdsP_,
                p,
                q
            );
        }
    }

    if (ionisationReactionQ)
    {
        ++nTotIonisationReactionsQ_;
        ++nIonisationReactionsQPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseQ
            (
                heatOfReactionIonQ_,
                ionisationProductsIdsQ_,
                p,
                q
            );
        }
    }

    // Perform exchange reaction
    if (exchangeReaction)
    {
        ++nTotExchangeReactions_;
        ++nExchangeReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;
            
            exchangePQ
            (
                heatOfReactionExchJoules,
                exchangeProductIds_,
                p,
                q
            );
        }
    }

    if (chargeExchangeReaction)
    {
        relax_ = false;

        nTotChargeExchangeReactions_++;
        nChargeExchangeReactionsPerTimeStep_++;

        if (allowSplitting_)
        {
            relax_ = false;
            
            scalar heatOfReactionChargeExchJoules =
                heatOfReactionChargeExchange_*physicoChemical::k.value();
            
            chargeExchangePQ
            (
                heatOfReactionChargeExchJoules,
                chargeExchangeProductIds_,
                p,
                q
            );
        }
    }
}


bool Foam::dissociationIonisationExchange::outputResults
(
    const label counterIndex
)
{
    label nTotExchangeReactions = nTotExchangeReactions_;
    label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
    label nTotDissociationReactions = nTotDissociationReactions_;
    label nTotIonisationReactionsP = nTotIonisationReactionsP_;
    label nTotIonisationReactionsQ = nTotIonisationReactionsQ_;
    
    label nDissociationReactionsPerTimeStep =
            nDissociationReactionsPerTimeStep_;
    label nIonisationReactionsPPerTimeStep =
        nIonisationReactionsPPerTimeStep_;
    label nIonisationReactionsQPerTimeStep =
        nIonisationReactionsQPerTimeStep_;
    label nExchangeReactionsPerTimeStep =
        nExchangeReactionsPerTimeStep_;
    label nChargeExchangeReactionsPerTimeStep =
        nChargeExchangeReactionsPerTimeStep_;
    
    if (Pstream::parRun())
    {
        reduce(nTotExchangeReactions, sumOp<label>());
        reduce(nTotChargeExchangeReactions, sumOp<label>());
        reduce(nTotDissociationReactions, sumOp<label>());
        reduce(nTotIonisationReactionsP, sumOp<label>());
        reduce(nTotIonisationReactionsQ, sumOp<label>());
        
        reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
        reduce(nIonisationReactionsPPerTimeStep, sumOp<label>());
        reduce(nIonisationReactionsQPerTimeStep, sumOp<label>());
        reduce(nExchangeReactionsPerTimeStep, sumOp<label>());
        reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
    }
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA =
        cloud_.typeIdList()[exchangeProductIds_[0]];
    const word& productB =
        cloud_.typeIdList()[exchangeProductIds_[1]];

    word productC;
    word productD;

    word productE;
    word productF;

    word productG;
    word productH;

    word productI;
    word productJ;
    
    if (!chargedMolecule_)
    {
        productC =
            cloud_.typeIdList()[dissociationProductIds_[0]];
        productD =
            cloud_.typeIdList()[dissociationProductIds_[1]];

        productE =
            cloud_.typeIdList()[ionisationProductsIdsP_[0]];
        productF =
            cloud_.typeIdList()[ionisationProductsIdsP_[1]];
    }

    if (!chargedAtom_)
    {
        productG =
            cloud_.typeIdList()[ionisationProductsIdsQ_[0]];
        productH =
            cloud_.typeIdList()[ionisationProductsIdsQ_[1]];
    }

    if (chargeExchange_)
    {
        productI =
            cloud_.typeIdList()[chargeExchangeProductIds_[0]];
        productJ =
            cloud_.typeIdList()[chargeExchangeProductIds_[1]];
    }
    
    bool write = uspReaction::outputResults(counterIndex);

    if (write)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols (2, 0);
        scalar volume = volume_;

        forAll(cellOccupancy, c)
        {
            const List<uspParcel*>& parcelsInCell = cellOccupancy[c];

            for (uspParcel* p : parcelsInCell)
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
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;
            scalar reactionRate5 = 0.0;

            reactionRate1 =
                (
                    nTotExchangeReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Exchange reaction "
                << reactantA << " + " << reactantB << " --> "
                << productA << " + " << productB
                << ", reaction rate = " << reactionRate1
                << endl;


            if (!chargedMolecule_)
            {
                reactionRate2 =
                    (
                        nTotDissociationReactions
                       *cloud_.nParticle()
                    )
                   /(
                        counterIndex*deltaT*numberDensities_[0]
                       *numberDensities_[1]*volume
                    );

                Info<< "Type II dissociation reaction "
                    << reactantA << " + " << reactantB << " --> "
                    << productC << " + " << productD
                    << " + " << reactantB
                    << ", reaction rate = " << reactionRate2
                    << endl;


                reactionRate3 =
                    (
                        nTotIonisationReactionsP
                       *cloud_.nParticle()
                    )
                   /(
                        counterIndex*deltaT*numberDensities_[0]
                       *numberDensities_[1]*volume
                    );

                Info<< "Ionisation reaction "
                    << reactantA << " + " << reactantB << " --> "
                    << productE << " + " << productF
                    << " + " << reactantB
                    << ", reaction rate = " << reactionRate3
                    << endl;
            }

            if (!chargedAtom_)
            {
                reactionRate4 =
                    (
                        nTotIonisationReactionsQ
                       *cloud_.nParticle()
                    )
                   /(
                        counterIndex*deltaT*numberDensities_[0]
                       *numberDensities_[1]*volume
                    );

                Info<< "Ionisation reaction "
                    << reactantA << " + " << reactantB << " --> "
                    << reactantA << " + "
                    << productG << " + " << productH
                    << ", reaction rate = " << reactionRate4
                    << endl;

            }

            if (chargeExchange_)
            {
                reactionRate5 =
                    (
                        nTotChargeExchangeReactions
                       *cloud_.nParticle()
                    )
                   /(
                        counterIndex*deltaT*numberDensities_[0]
                       *numberDensities_[1]*volume
                    );

                Info<< "Charge exchange reaction "
                    << reactantA << " + " << reactantB << " --> "
                    << productI << " + " << productJ
                    << ", reaction rate = " << reactionRate5
                    << endl;
            }
        }
    }
    else
    {
        if (nTotExchangeReactions > VSMALL)
        {
            Info<< "Exchange reaction "
                << reactantA << " + " << reactantB << " --> "
                << productA << " + " << productB
                << " is active, nReactions this time step = "
                << nExchangeReactionsPerTimeStep << endl;
        }

        if (nTotDissociationReactions > VSMALL)
        {
            Info<< "Type II dissociation reaction "
                << reactantA << " + " << reactantB << " --> "
                << productC << " + " << productD
                << " + " << reactantB
                << " is active, nReactions this time step = "
                << nDissociationReactionsPerTimeStep << endl;
        }

        if (nTotIonisationReactionsP > VSMALL)
        {
            Info<< "Ionisation reaction "
                << reactantA << " + " << reactantB << " --> "
                << productE << " + " << productF
                << " + " << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPPerTimeStep << endl;
        }

        if (nTotIonisationReactionsQ > VSMALL)
        {
            Info<< "Ionisation reaction "
                << reactantA << " + " << reactantB << " --> "
                << reactantA << " + "
                << productG << " + " << productH
                << " is active, nReactions this time step = "
                << nIonisationReactionsQPerTimeStep << endl;
        }

        if (nTotChargeExchangeReactions > VSMALL)
        {
            Info<< "Charge exchange reaction "
                << reactantA << " + " << reactantB << " --> "
                << productI << " + " << productJ
                << " is active, nReactions this time step = "
                << nChargeExchangeReactionsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPPerTimeStep_ = 0.0;
    nIonisationReactionsQPerTimeStep_ = 0.0;
    nExchangeReactionsPerTimeStep_ = 0.0;
    nChargeExchangeReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
