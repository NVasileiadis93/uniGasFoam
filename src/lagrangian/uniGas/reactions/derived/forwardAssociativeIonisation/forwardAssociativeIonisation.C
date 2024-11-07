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

#include "forwardAssociativeIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(forwardAssociativeIonisation, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    forwardAssociativeIonisation,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forwardAssociativeIonisation::forwardAssociativeIonisation
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    intermediateId_(),
    associativeIonisationProductIds_(),
    ionisationProductIds_(),
    heatOfReactionRecombination_(),
    heatOfReactionIntermediateIonisation_(),
    heatOfReactionIonisation_(),
    nTotalIonisationReactions_(0),
    nTotalAssociativeIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    nAssociativeIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forwardAssociativeIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::forwardAssociativeIonisation::setProperties()
{
    // reading in reactants

    forAll(reactantIds_, r)
    {
        // check that reactants are 'atoms'

        const label rDofReactant =
            cloud_.constProps(reactantIds_[r]).rotationalDoF();

        if (rDofReactant > VSMALL)
        {
            FatalErrorInFunction
                << "Reactant must be an atom (not a molecule): "
                << reactants_[r]
                << exit(FatalError);
        }
    }

    // reading in associative ionisation products

    const wordPair associativeIonisationProductMolecules
    (
        propsDict_.lookup("productsOfAssociativeIonisation")
    );

    associativeIonisationProductIds_ = labelPair(-1, -1);

    forAll(associativeIonisationProductIds_, i)
    {
        associativeIonisationProductIds_[i] =
            cloud_.typeIdList().find(associativeIonisationProductMolecules[i]);

        // check that products belong to the typeIdList
        if (associativeIonisationProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: "
                << associativeIonisationProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are a 'molecule' and an 'electron'

    labelPair assIonProds = associativeIonisationProductIds_;

    label rDofProd1 = cloud_.constProps(assIonProds[0]).rotationalDoF();

    if (rDofProd1 < 1)
    {
        FatalErrorInFunction
            << "First product must be a molecule: "
            << associativeIonisationProductMolecules
            << exit(FatalError);
    }

    const label charge =
        cloud_.constProps(associativeIonisationProductIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second product must be an electron: "
            << associativeIonisationProductMolecules
            << exit(FatalError);
    }

    // reading in ionisation of P products

    const wordPair ionisationProductMolecules
    (
        propsDict_.lookup("productsOfIonisation")
    );

    ionisationProductIds_ = labelPair(-1, -1);

    forAll(ionisationProductIds_, i)
    {
        ionisationProductIds_[i] =
            cloud_.typeIdList().find(ionisationProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << ionisationProductMolecules[i]
                << exit(FatalError);
        }
    }


    // check that products are an 'atom' and an 'electron'

    rDofProd1 = cloud_.constProps(ionisationProductIds_[0]).rotationalDoF();

    label charge2 = cloud_.constProps(ionisationProductIds_[0]).charge();

    if (rDofProd1 > VSMALL || charge2 == -1)
    {
        FatalErrorInFunction
            << "First product must be a charged atom: "
            << ionisationProductMolecules
            << exit(FatalError);
    }

    label charge3 = cloud_.constProps(ionisationProductIds_[1]).charge();

    if (charge3 != -1)
    {
        FatalErrorInFunction
            << "Second product must be an electron: "
            << ionisationProductMolecules
            << exit(FatalError);
    }

    // reading in intermediate molecule

    word intermediateMolecule(propsDict_.get<word>("intermediateMolecule"));

    intermediateId_ = cloud_.typeIdList().find(intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/uniGasProperties)
    if (intermediateId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find type id: " << intermediateMolecule << nl
            << exit(FatalError);
    }

    // check that the intermediate is a 'molecule'

    const label rDof = cloud_.constProps(intermediateId_).rotationalDoF();

    if (rDof < 1)
    {
        FatalErrorInFunction
            << "The intermediate specie must be a molecule (not an atom): "
            << intermediateMolecule
            << exit(FatalError);
    }

    heatOfReactionRecombination_ =
        propsDict_.get<scalar>("heatOfReactionRecombination");
    heatOfReactionIntermediateIonisation_ =
        propsDict_.get<scalar>("heatOfReactionIntermediateIonisation");
    heatOfReactionIonisation_ =
        propsDict_.get<scalar>("heatOfReactionIonisation");
}


bool Foam::forwardAssociativeIonisation::tryReactMolecules
(
    const label typeIdP,
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


void Foam::forwardAssociativeIonisation::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    const label typeIdP = p.typeId();

    if (reactantIds_[0] == reactantIds_[1] && typeIdP == reactantIds_[0])
    {
        react
        (
            p,
            q          
        );
    }
}

void Foam::forwardAssociativeIonisation::react
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    relax_ = true;

    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(2, 0.0);
    const auto& constPropInter = cloud_.constProps(intermediateId_);

    scalar omegaIntermediate = constPropInter.omega();
    label rotationalDofIntermediate = constPropInter.rotationalDoF();
    scalar ChiBIntermediate = 2.5 - omegaIntermediate;
    scalar thetaVIntermediate = constPropInter.thetaV()[0];
    scalar thetaDIntermediate = constPropInter.thetaD()[0];
    scalar ZrefIntermediate = constPropInter.Zref()[0];
    scalar refTempZvIntermediate = constPropInter.TrefZv()[0];
    label ELevelIntermediate = 0;
    const List<scalar>& EElistIntermediate =
        constPropInter.electronicEnergyList();
    const List<label>& gListIntermediate = constPropInter.degeneracyList();
    label jMaxIntermediate = constPropInter.nElectronicLevels();

    bool ionisationReaction = false;
    bool associativeIonisation = false;

    // 2 reactions possible

    // 1. Ionisation of P
    // 2. Forward associative ionisation

    scalar Ec = 0.0;

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of P is possible
    Ec = translationalEnergy(p, q) + EEle(p);

    if ((Ec - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[0] = 1.0;
    }

    // collision energy is the translational energy of the two atoms,
    // plus their electronic energies

    Ec = translationalEnergy(p, q) + EEle(p);
    scalar EcOrig = Ec;

    label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

    label vibLevelIntermediate = -1;

    if (iMax > SMALL)
    {
        vibLevelIntermediate = cloud_.postCollisionVibrationalEnergyLevel
        (
            true,
            0.0,
            iMax,
            thetaVIntermediate,
            thetaDIntermediate,
            refTempZvIntermediate,
            omegaIntermediate,
            ZrefIntermediate,
            Ec
            );
    }

    if (vibLevelIntermediate == 0)
    {
        //'Form' the intermediate molecule and test it for ionisation
        const scalar heatOfReactionRecombinationJoules =
            heatOfReactionRecombination_*physicoChemical::k.value();

        Ec = EcOrig + heatOfReactionRecombinationJoules;

        label postCollisionELevel =
            cloud_.postCollisionElectronicEnergyLevel
            (
                Ec,
                jMaxIntermediate,
                omegaIntermediate,
                EElistIntermediate,
                gListIntermediate
            );

        ELevelIntermediate = postCollisionELevel;

        // relative translational energy after electronic exchange
        Ec -= EElistIntermediate[ELevelIntermediate];

        iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

        if (iMax > SMALL)
        {
            label postCollisionVibLevel =
                cloud_.postCollisionVibrationalEnergyLevel
                (
                    true,
                    0.0,
                    iMax,
                    thetaVIntermediate,
                    thetaDIntermediate,
                    refTempZvIntermediate,
                    omegaIntermediate,
                    ZrefIntermediate,
                    Ec
                    );

            Ec -=
                postCollisionVibLevel
                *thetaVIntermediate*physicoChemical::k.value();
        }

        scalar ERot = 0.0;

        scalar energyRatio =
            cloud_.postCollisionRotationalEnergy
            (
                rotationalDofIntermediate,
                ChiBIntermediate
            );

        ERot = energyRatio*Ec;

        Ec -= ERot;

        // redistribution finished, test it for ionisation

        scalar EcPPIon = 0.0;
        scalar ionisationEnergy =
            cloud_.constProps(intermediateId_).ionisationTemperature()
            *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        EcPPIon = Ec + EElistIntermediate[ELevelIntermediate];

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            // Associative ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        scalarList probNorm
        (
            reactionProbabilities.size(),
            0.0
        );

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
                        ionisationReaction = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Associative ionisation reaction is to occur
                        associativeIonisation = true;
                        break;
                    }
                }
            }
        }
    }

    if (ionisationReaction)
    {
        ++nTotalIonisationReactions_;
        ++nIonisationReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonisation_,
                ionisationProductIds_,
                p,
                q
            );
        }
    }

    if (associativeIonisation)
    {
        nTotalAssociativeIonisationReactions_++;
        nAssociativeIonisationReactionsPerTimeStep_++;

        if (allowSplitting_)
        {
            relax_ = false;

            this->associativeIonisation
            (
                heatOfReactionIntermediateIonisation_,
                heatOfReactionRecombination_,
                associativeIonisationProductIds_,
                p,
                q
            );
        }
    }
}


bool Foam::forwardAssociativeIonisation::outputResults(const label counterIndex)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA =
        cloud_.typeIdList()[associativeIonisationProductIds_[0]];
    const word& productB =
        cloud_.typeIdList()[associativeIonisationProductIds_[1]];

    const word& productC = cloud_.typeIdList()[ionisationProductIds_[0]];
    const word& productD = cloud_.typeIdList()[ionisationProductIds_[1]];

    if (write)
    {
        // measure density

        const auto& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols (2, 0);
        scalar volume = volume_;
        label nTotalAssociativeIonisationReactions =
            nTotalAssociativeIonisationReactions_;
        label nTotalIonisationReactions = nTotalIonisationReactions_;

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
        }

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationReactions, sumOp<label>());

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

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateAssociativeIonisation = 0.0;
            scalar reactionRateIonisation = 0.0;


            reactionRateAssociativeIonisation =
                (
                    nTotalAssociativeIonisationReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Associative ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB
                << ", reaction rate = " << reactionRateAssociativeIonisation
                << endl;

            reactionRateIonisation =
                (
                    nTotalIonisationReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productC << " + " << productD << " + " << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotalAssociativeIonisationReactions =
            nTotalAssociativeIonisationReactions_;
        label nTotalIonisationReactions = nTotalIonisationReactions_;

        label nAssociativeIonisationReactionsPerTimeStep =
            nAssociativeIonisationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep =
            nIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationReactions, sumOp<label>());

            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotalAssociativeIonisationReactions > VSMALL)
        {
            Info<< "Associative ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB
                << " is active, nReactions this time step = "
                << nAssociativeIonisationReactionsPerTimeStep << endl;
        }

        if (nTotalIonisationReactions > VSMALL)
        {
            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nIonisationReactionsPerTimeStep << endl;
        }
    }

    nAssociativeIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
