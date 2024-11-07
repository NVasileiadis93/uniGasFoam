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

#include "atomElectronIonisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(atomElectronIonisation, 0);

addToRunTimeSelectionTable
(
    uniGasReaction,
    atomElectronIonisation,
    dictionary
);
}


void Foam::atomElectronIonisation::setProperties()
{
    if (reactants_[0] == reactants_[1])
    {
        FatalErrorInFunction
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    // check that reactant one is an 'atom'

    const label rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDoF();

    const label vDof1 = cloud_.constProps(reactantIds_[0]).vibrationalDoF();

    if (rDof1 > VSMALL)
    {
        FatalErrorInFunction
            << "First reactant must be an atom "
            << "(not a molecule or an electron): " << reactants_[0]
            << exit(FatalError);
    }

    if (vDof1 > VSMALL)
    {
         FatalErrorInFunction
            << "Reactions are currently only implemented "
            << "for monatomic and diatomic species."
            << " This is a polyatomic:" << reactants_[0]
            << exit(FatalError);
    }

    const label charge1 = cloud_.constProps(reactantIds_[0]).charge();

    if (charge1 == -1)
    {
        FatalErrorInFunction
            << "First reactant must not be an electron "
            << reactants_[1]
            << exit(FatalError);
    }

    // check that reactant two is an 'electron'

    const label charge2 = cloud_.constProps(reactantIds_[1]).charge();

    if (charge2 != -1)
    {
        FatalErrorInFunction
            << "Second reactant must be an electron, not "
            << reactants_[1]
            << exit(FatalError);
    }

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

        // check that product one is a 'atom'

        const label rDof2 =
            cloud_.constProps(productIdsIon_[0]).rotationalDoF();

        const label vDof2 =
            cloud_.constProps(productIdsIon_[0]).vibrationalDoF();

        if (rDof2 > 1)
        {
            FatalErrorInFunction
                << "First product must be an atom (not an atom): "
                << productsOfIonisedAtom[0]
                << exit(FatalError);
        }

        if (vDof2 > VSMALL)
        {
            FatalErrorInFunction
                << "Reactions are currently only implemented "
                << "for monatomic and diatomic species."
                << " This is a polyatomic:" << productsOfIonisedAtom[0]
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

    heatOfReactionIon_ = propsDict_.get<scalar>("heatOfReactionIonisation");
}


void Foam::atomElectronIonisation::react
(
    const labelPair& productIDs,
    uniGasParcel& p,
    uniGasParcel& q
)
{
    relax_ = true;

    scalar ionisationEnergy =
        cloud_.constProps(p.typeId()).ionisationTemperature()
        *physicoChemical::k.value();

    scalar EcPPIon = translationalEnergy(p, q) + EEle(p);

    if ((EcPPIon - ionisationEnergy) > VSMALL)
    {
        ++nIonisationReactions_;
        ++nReactionsPerTimeStep_;

        if (allowSplitting_)
        {
            relax_ = false;
            ioniseP(heatOfReactionIon_, productIDs, p, q);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atomElectronIonisation::atomElectronIonisation
(
    const Time& t,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasReaction(t, cloud, dict),
    productIdsIon_(),
    heatOfReactionIon_(),
    nIonisationReactions_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atomElectronIonisation::initialConfiguration()
{
    setCommonReactionProperties();
    setProperties();
}


bool Foam::atomElectronIonisation::tryReactMolecules
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


void Foam::atomElectronIonisation::reaction
(
    uniGasParcel& p,
    uniGasParcel& q
)
{
    if (p.typeId() == reactantIds_[0] && q.typeId() == reactantIds_[1])
    {
        react(productIdsIon_, p, q);
    }
    else if (p.typeId() == reactantIds_[1] && q.typeId() == reactantIds_[0])
    {
        react(productIdsIon_, q, p);
    }
}


bool Foam::atomElectronIonisation::outputResults(const label counterIndex)
{
    bool write = uniGasReaction::outputResults(counterIndex);
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
    
    const word& productA = cloud_.typeIdList()[productIdsIon_[0]];
    const word& productB = cloud_.typeIdList()[productIdsIon_[1]];

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
                    ++mols[id];
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        label nTotReactionsIonisation = nIonisationReactions_;

        // Parallel communication
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;

            reactionRateIonisation =
            (
                nTotReactionsIonisation
               *cloud_.nParticle()
            )
           /(
                counterIndex*deltaT*numberDensities_[0]
               *numberDensities_[1]*volume
            );

            Info<< "Electron ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotReactionsIonisation = nIonisationReactions_;
        label nReactionsPerTimeStep = nReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotReactionsIonisation > VSMALL)
        {
            Info<< "Electron ionisation reaction "
                <<  reactantA << " + " << reactantB
                <<  " --> "
                << productA << " + " << productB << " + "
                << reactantB
                << " is active, nReactions this time step = "
                << nReactionsPerTimeStep << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;

    return write;
}


// ************************************************************************* //
