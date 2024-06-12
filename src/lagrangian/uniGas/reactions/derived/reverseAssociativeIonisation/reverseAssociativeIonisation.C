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

#include "reverseAssociativeIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(reverseAssociativeIonisation, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    reverseAssociativeIonisation,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reverseAssociativeIonisation::reverseAssociativeIonisation
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    intermediateId_(),
    productIds_(),
    heatOfReactionDissociation_(),
    heatOfReactionRecombination_(),
    nReactions_(),
    nReactionsPerTimeStep_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reverseAssociativeIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


void Foam::reverseAssociativeIonisation::setProperties()
{
    // check that first reactant is a charged molecule

    const label rDofReactant1 =
        cloud_.constProps(reactantIds_[0]).rotationalDoF();
    const label chargeReactant1 =
        cloud_.constProps(reactantIds_[0]).charge();

    if (rDofReactant1 < VSMALL)
    {
        FatalErrorInFunction
            << "First reactant must be an ionised molecule: "
            << reactants_[0]
            << exit(FatalError);
    }

    if (chargeReactant1 != 1)
    {
        FatalErrorInFunction
            << "First reactant must be an ionised molecule: "
            << reactants_[0]
            << exit(FatalError);
    }

    const label vDofReactant1 =
        cloud_.constProps(reactantIds_[0]).vibrationalDoF();

    if (vDofReactant1 > 1)
    {
        FatalErrorInFunction
            << "Reactions are currently only implemented "
            << "for monatomic and diatomic species"
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    // check that second reactant is an 'electron'

    const label charge = cloud_.constProps(reactantIds_[1]).charge();

    if (charge != -1)
    {
        FatalErrorInFunction
            << "Second reactant must be an electron: " << reactants_[0]
            << exit(FatalError);
    }

    // reading in product

    const wordPair productMolecules
    (
        propsDict_.lookup("productsOfAssociativeIonisation")
    );

    productIds_ = labelPair(-1, -1);

    forAll(productIds_, i)
    {
        productIds_[i] = cloud_.typeIdList().find(productMolecules[i]);

        // check that products belong to the typeIdList
        if (productIds_[i] == -1)
        {
            FatalErrorInFunction
                << "Cannot find type id: " << productMolecules[i]
                << exit(FatalError);
        }

        // check that products are 'atoms'

        label rDof = cloud_.constProps(productIds_[i]).rotationalDoF();

        if (rDof > 1)
        {
            FatalErrorInFunction
                << "Product must be an atom: " << productMolecules[i]
                << exit(FatalError);
        }
    }


    // reading in intermediate molecule

    const word intermediateMolecule
    (
        propsDict_.get<word>("intermediateMolecule")
    );

    intermediateId_ = cloud_.typeIdList().find(intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/uniGasProperties)
    if (intermediateId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find type id: " << intermediateMolecule
            << exit(FatalError);
    }

    // check that the intermediate is a 'molecule'

    label rDof = cloud_.constProps(intermediateId_).rotationalDoF();

    if (rDof < 1)
    {
        FatalErrorInFunction
            << "The intermediate specie must be a molecule "
            << "(not an atom): " << intermediateMolecule
            << exit(FatalError);
    }

    heatOfReactionDissociation_ =
        propsDict_.get<scalar>("heatOfReactionDissociation");
    heatOfReactionRecombination_ =
        propsDict_.get<scalar>("heatOfReactionRecombination");
}


bool Foam::reverseAssociativeIonisation::tryReactMolecules
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


void Foam::reverseAssociativeIonisation::reaction
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

void Foam::reverseAssociativeIonisation::react
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    // Perform a 'recombination' to test the intermediate molecule

    relax_ = true;

    const auto& constPropsInter = cloud_.constProps(intermediateId_);

    scalar omegaIntermediate = constPropsInter.omega();
    label rotationalDofIntermediate = constPropsInter.rotationalDoF();
    scalar ChiBIntermediate = 2.5 - omegaIntermediate;
    scalar thetaVIntermediate = constPropsInter.thetaV()[0];
    scalar thetaDIntermediate = constPropsInter.thetaD()[0];
    scalar ZrefIntermediate = constPropsInter.Zref()[0];
    scalar refTempZvIntermediate = constPropsInter.TrefZv()[0];
    label ELevelIntermediate = -1;
    const List<scalar>& EElistIntermediate =
        constPropsInter.electronicEnergyList();
    const List<label>& gListIntermediate = constPropsInter.degeneracyList();
    label jMaxIntermediate = constPropsInter.nElectronicLevels();

    // collision energy is the translational energy of the two atoms,
    // plus their electronic energies

    scalar Ec = translationalEnergy(p, q) + EEle(p);
    scalar EcOrig = Ec;

    label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

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

    if (ELevelIntermediate == 0)
    {
        // 'Form' the intermediate molecule and test it for ionisation
        const scalar heatOfReactionRecombinationJoules =
            heatOfReactionRecombination_*physicoChemical::k.value();

        Ec = EcOrig + heatOfReactionRecombinationJoules;

        ELevelIntermediate =
            cloud_.postCollisionElectronicEnergyLevel
            (
                Ec,
                jMaxIntermediate,
                omegaIntermediate,
                EElistIntermediate,
                gListIntermediate
            );

        iMax = (Ec/(physicoChemical::k.value()*thetaVIntermediate));

        label eVibLevel = -1;

        if (iMax > SMALL)
        {
            eVibLevel =
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

            Ec -= eVibLevel*physicoChemical::k.value()*thetaVIntermediate;
        }

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

        // redistribution finished, test it for dissociation

        scalar EcDiss = 0.0;
        label iMax = 0;
        label id =
            cloud_.constProps(intermediateId_).charDissQuantumLevel()[0];

        // calculate if a dissociation is possible
        scalar EVib =
            eVibLevel*physicoChemical::k.value()*thetaVIntermediate;

        EcDiss = Ec + EVib;
        iMax = EcDiss/(physicoChemical::k.value()*thetaVIntermediate);

        if ((iMax - id) > VSMALL)
        {
            ++nReactions_;
            ++nReactionsPerTimeStep_;

            if (allowSplitting_)
            {
                relax_ = false;

                associativeIonisation
                (
                    heatOfReactionRecombination_,
                    heatOfReactionDissociation_,
                    productIds_,
                    p,
                    q
                );
            }
        }
    }
}


bool Foam::reverseAssociativeIonisation::outputResults(const label counterIndex)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIds_[0]];
    const word& productB = cloud_.typeIdList()[productIds_[1]];

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
        label nTotReactions = nReactions_;

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
            reduce(nTotReactions, sumOp<label>());
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
            scalar reactionRate = 0.0;

            reactionRate =
                (nTotReactions*cloud_.nParticle())
               /(
                    counterIndex*deltaT*numberDensities_[0]
                    *numberDensities_[1]*volume
                );

            Info<< "Associative ionisation reaction "
                << reactantA << " + " << reactantB
                << " --> "
                << productA << " + " << productB
                << ", reaction rate = " << reactionRate
                << endl;
        }
    }
    else
    {
        label nTotReactions = nReactions_;
        label nReactionsPerTimeStep = nReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactions, sumOp<label>());
            reduce(nReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotReactions > VSMALL)
        {
            Info<< "Associative ionisation reaction "
                << reactantA << " + " << reactantB
                << " --> "
                << productA << " + " << productB
                << " is active, nReactions this time step = "
                << nReactionsPerTimeStep << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
