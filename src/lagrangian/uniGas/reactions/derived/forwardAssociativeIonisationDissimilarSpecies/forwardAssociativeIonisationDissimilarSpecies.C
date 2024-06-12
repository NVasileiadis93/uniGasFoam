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

#include "forwardAssociativeIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(forwardAssociativeIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    forwardAssociativeIonisationDissimilarSpecies,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forwardAssociativeIonisationDissimilarSpecies::
forwardAssociativeIonisationDissimilarSpecies
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    intermediateId_(),
    associativeIonisationProductIds_(),
    ionisationPProductIds_(),
    ionisationQProductIds_(),
    heatOfReactionRecombination_(),
    heatOfReactionIntermediateIonisation_(),
    heatOfReactionIonisationP_(),
    heatOfReactionIonisationQ_(),
    nTotalIonisationPReactions_(0),
    nTotalIonisationQReactions_(0),
    nTotalAssociativeIonisationReactions_(0),
    nIonisationPReactionsPerTimeStep_(0),
    nIonisationQReactionsPerTimeStep_(0),
    nAssociativeIonisationReactionsPerTimeStep_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forwardAssociativeIonisationDissimilarSpecies::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::forwardAssociativeIonisationDissimilarSpecies::setProperties()
{
    forAll(reactantIds_, r)
    {
        // check that reactants are 'atoms'

        label rDofReactant =
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
                << associativeIonisationProductMolecules[i]
                << exit(FatalError);
        }
    }


    // check that products are a 'molecule' and an 'electron'

    label id0 = associativeIonisationProductIds_[0];

    label rDofProd1 = cloud_.constProps(id0).rotationalDoF();

    if (rDofProd1 < 1)
    {
        FatalErrorInFunction
            << "First product must be a molecule: "
            << associativeIonisationProductMolecules
            << exit(FatalError);
    }

    label id1 = associativeIonisationProductIds_[1];

    label charge = cloud_.constProps(id1).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second product must be an electron: "
            << associativeIonisationProductMolecules
            << exit(FatalError);
    }

    // reading in ionisation of P products

    const wordPair ionisationPProductMolecules
    (
        propsDict_.lookup("productsOfIonisationP")
    );

    ionisationPProductIds_ = labelPair(-1, -1);

    forAll(ionisationPProductIds_, i)
    {
        ionisationPProductIds_[i] =
            cloud_.typeIdList().find(ionisationPProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationPProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: "
                << ionisationPProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are an 'atom' and an 'electron'

    rDofProd1 = cloud_.constProps(ionisationPProductIds_[0]).rotationalDoF();

    if (rDofProd1 > 0)
    {
        FatalErrorInFunction
            << "First product must be a charged atom: "
            << ionisationPProductMolecules
            << exit(FatalError);
    }

    const label charge2 = cloud_.constProps(ionisationPProductIds_[1]).charge();

    if (charge2 != -1)
    {
        FatalErrorInFunction
            << "Second product must be an electron: "
            << ionisationPProductMolecules
            << exit(FatalError);
    }

    // reading in ionisation of Q products

    const wordPair ionisationQProductMolecules
    (
        propsDict_.lookup("productsOfIonisationQ")
    );

    ionisationQProductIds_ = labelPair(-1, -1);

    forAll(ionisationQProductIds_, i)
    {
        ionisationQProductIds_[i] =
            cloud_.typeIdList().find(ionisationQProductMolecules[i]);

        // check that products belong to the typeIdList
        if (ionisationQProductIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: "
                << ionisationQProductMolecules[i] << nl
                << exit(FatalError);
        }
    }


    // check that products are an 'atom' and an 'electron'

    rDofProd1 = cloud_.constProps(ionisationQProductIds_[0]).rotationalDoF();

    if (rDofProd1 > 0)
    {
        FatalErrorInFunction
            << "First product must be a charged atom: "
            << ionisationQProductMolecules
            << exit(FatalError);
    }

    const label charge3 = cloud_.constProps(ionisationQProductIds_[1]).charge();

    if (charge3 != -1)
    {
        FatalErrorInFunction
            << "Second product must be an electron: "
            << ionisationQProductMolecules
            << exit(FatalError);
    }

    // reading in intermediate molecule

    const word intermediateMolecule =
        propsDict_.get<word>("intermediateMolecule");

    intermediateId_ = cloud_.typeIdList().find(intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/uniGasProperties)
    if (intermediateId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find type id: "
            << intermediateMolecule
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
    heatOfReactionIonisationP_ =
        propsDict_.get<scalar>("heatOfReactionIonisationP");
    heatOfReactionIonisationQ_ =
        propsDict_.get<scalar>("heatOfReactionIonisationQ");
}


bool Foam::forwardAssociativeIonisationDissimilarSpecies::tryReactMolecules
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


void Foam::forwardAssociativeIonisationDissimilarSpecies::reaction
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
    else if (typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        react
        (
            q,
            p          
        );
    }
}

void Foam::forwardAssociativeIonisationDissimilarSpecies::react
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    relax_ = true;

    scalar totalReactionProbability = 0.0;
    scalarList reactionProbabilities(3, 0.0);

    const auto& constPropsInter = cloud_.constProps(intermediateId_);
    const scalar omegaIntermediate = constPropsInter.omega();
    const label rotationalDofIntermediate = constPropsInter.rotationalDoF();
    const scalar ChiBIntermediate = 2.5 - omegaIntermediate;
    const scalar thetaVIntermediate = constPropsInter.thetaV()[0];
    const scalar thetaDIntermediate = constPropsInter.thetaD()[0];
    const scalar ZrefIntermediate = constPropsInter.Zref()[0];
    const scalar refTempZvIntermediate = constPropsInter.TrefZv()[0];
    label ELevelIntermediate = 0;
    const List<scalar>& EElistIntermediate =
        constPropsInter.electronicEnergyList();
    const List<label>& gListIntermediate =
        constPropsInter.degeneracyList();
    const label jMaxIntermediate = constPropsInter.nElectronicLevels();


    bool ionisationReactionP = false;
    bool ionisationReactionQ = false;
    bool associativeIonisation = false;

    // 3 reactions possible

    // 1. Ionisation of P
    // 2. Ionisation of Q
    // 3. Forward associative ionisation

    scalar Ec = 0.0;

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of species P is possible
    Ec = translationalEnergy(p, q) + EEle(p);

    if ((Ec - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[0] = 1.0;
    }

    ionisationEnergy =
        cloud_.constProps(q.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    // calculate if an ionisation of species Q is possible
    Ec = translationalEnergy(p, q) + EEle(q);

    if ((Ec - ionisationEnergy) > VSMALL)
    {
        // Ionisation can occur
        totalReactionProbability += 1.0;
        reactionProbabilities[1] = 1.0;
    }

    Ec = translationalEnergy(p, q) + EEle(p);
    scalar EcOrig = Ec;

    label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

    label vibLevelIntermediate = -1;

    if (iMax > SMALL)
    {
        vibLevelIntermediate =
            cloud_.postCollisionVibrationalEnergyLevel
            (
                false,
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
                postCollisionVibLevel*thetaVIntermediate
                *physicoChemical::k.value();
        }

        const label postCollisionELevel =
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

        const scalar ionisationEnergy =
            cloud_.constProps(intermediateId_).ionisationTemperature()
            *physicoChemical::k.value();

        // calculate if an ionisation of species P is possible
        const scalar EcPPIon = Ec + EElistIntermediate[ELevelIntermediate];

        if ((EcPPIon - ionisationEnergy) > VSMALL)
        {
            // Associative ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }
    }

    // Decide if a reaction is to occur

    if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
    {
        // A chemical reaction is to occur, choose which one

        const scalarList probNorm(reactionProbabilities/totalReactionProbability);

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
                        // Ionisation is to occur
                        ionisationReactionP = true;
                        break;
                    }
                    if (i == 1)
                    {
                        // Ionisation reaction is to occur
                        ionisationReactionQ = true;
                        break;
                    }
                    if (i == 2)
                    {
                        // Associative ionisation is to occur
                        associativeIonisation = true;
                        break;
                    }
                }
            }
        }
    }

    if (ionisationReactionP)
    {
        ++nTotalIonisationPReactions_;
        ++nIonisationPReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseP
            (
                heatOfReactionIonisationP_,
                ionisationPProductIds_,
                p,
                q
            );
        }
    }

    if (ionisationReactionQ)
    {
        ++nTotalIonisationQReactions_;
        ++nIonisationQReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;

            ioniseQ
            (
                heatOfReactionIonisationQ_,
                ionisationQProductIds_,
                p,
                q
            );
        }
    }

    if (associativeIonisation)
    {
        ++nTotalAssociativeIonisationReactions_;
        ++nAssociativeIonisationReactionsPerTimeStep_;

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

bool Foam::forwardAssociativeIonisationDissimilarSpecies::outputResults
(
    const label counterIndex
)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA =
        cloud_.typeIdList()[associativeIonisationProductIds_[0]];
    const word& productB =
        cloud_.typeIdList()[associativeIonisationProductIds_[1]];

    const word& productC =
        cloud_.typeIdList()[ionisationPProductIds_[0]];
    const word& productD =
        cloud_.typeIdList()[ionisationPProductIds_[1]];

    const word& productE =
        cloud_.typeIdList()[ionisationQProductIds_[0]];
    const word& productF =
        cloud_.typeIdList()[ionisationQProductIds_[1]];

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
        label nTotalAssociativeIonisationReactions =
            nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;

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
        }

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());
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
            const scalar reactionRateAssociativeIonisation =
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

            const scalar reactionRateIonisationP =
                (
                    nTotalIonisationPReactions
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
                << ", reaction rate = " << reactionRateIonisationP
                << endl;

            const scalar reactionRateIonisationQ =
                (
                    nTotalIonisationQReactions
                   *cloud_.nParticle()
                )
               /(
                    counterIndex*deltaT*numberDensities_[0]
                   *numberDensities_[1]*volume
                );

            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << reactantA << " + " << productE << " + " << productF
                << ", reaction rate = " << reactionRateIonisationQ
                << endl;
        }
    }
    else
    {
        label nTotalAssociativeIonisationReactions =
            nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;

        label nAssociativeIonisationReactionsPerTimeStep =
            nAssociativeIonisationReactionsPerTimeStep_;
        label nIonisationPReactionsPerTimeStep =
            nIonisationPReactionsPerTimeStep_;
        label nIonisationQReactionsPerTimeStep =
            nIonisationQReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());

            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationPReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationQReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotalAssociativeIonisationReactions > VSMALL)
        {
            Info<< "Associative ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB
                << " is active, nReactions this time step = "
                << nAssociativeIonisationReactionsPerTimeStep
                << endl;
        }

        if (nTotalIonisationPReactions > VSMALL)
        {
            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nIonisationPReactionsPerTimeStep
                << endl;
        }

        if (nTotalIonisationQReactions > VSMALL)
        {
            Info<< "Ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << reactantA << " + " << productA
                << " + " << productB
                << " is active, nReactions this time step = "
                << nIonisationQReactionsPerTimeStep
                << endl;
        }
    }

    nAssociativeIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationPReactionsPerTimeStep_ = 0.0;
    nIonisationQReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
