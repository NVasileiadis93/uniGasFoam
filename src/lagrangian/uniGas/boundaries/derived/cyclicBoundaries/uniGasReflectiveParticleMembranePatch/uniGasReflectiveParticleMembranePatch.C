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

#include "uniGasReflectiveParticleMembranePatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasReflectiveParticleMembranePatch, 0);

addToRunTimeSelectionTable
(
    uniGasCyclicBoundary,
    uniGasReflectiveParticleMembranePatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasReflectiveParticleMembranePatch::uniGasReflectiveParticleMembranePatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasCyclicBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    p_(propsDict_.get<scalar>("reflectionProbability")),
    temperature_(propsDict_.get<scalar>("temperature")),
    velocity_(propsDict_.get<vector>("velocity")),
    nReflections_(0),
    nRejections_(0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasReflectiveParticleMembranePatch::calculateProperties()
{
    label nReflections = nReflections_;
    label nRejections = nRejections_;

    if (Pstream::parRun())
    {
        reduce(nReflections, sumOp<label>());
        reduce(nRejections, sumOp<label>());
    }

    if (nRejections > 0)
    {
        Info<< "    no Reflections: " << nReflections
            << ", no Rejections: " << nRejections
            << " ratio reflections/(reflections+rejections):"
            << scalar(nReflections)/scalar(nReflections+nRejections)
            << endl;
    }
}


void Foam::uniGasReflectiveParticleMembranePatch::initialConfiguration()
{}


void Foam::uniGasReflectiveParticleMembranePatch::controlMol
(
    uniGasParcel& p,
    uniGasParcel::trackingData& td
)
{
    const label faceI = p.face();

    vector nF = mesh_.faceAreas()[faceI];
    vector nw = p.normal();

    vector& U = p.U();

    label fA = coupledFacesA_.find(faceI);
    label fB = coupledFacesB_.find(faceI);

    nF /= mag(nF);
    nw /= mag(nw);

    scalar U_dot_nw = U & nw;
    vector Ut = U - U_dot_nw*nw;

    label typeId = p.typeId();

    const scalar T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    Random& rndGen = cloud_.rndGen();

    scalar d = nF & U;

    if(d > 0) // processor boundary
    {
        if(fA != -1)
        {
            scalar pRandom = rndGen.sample01<scalar>();

            if( pRandom <= p_ ) // reflect molecule
            {
                while (mag(Ut) < SMALL)
                {
                    // If the incident velocity is parallel to the face normal, no
                    // tangential direction can be chosen.  Add a perturbation to the
                    // incoming velocity and recalculate.

                    U = vector
                    (
                        U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                        U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                        U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
                    );

                    U_dot_nw = U & nw;

                    Ut = U - U_dot_nw*nw;
                }

                // Wall tangential unit vector
                vector tw1 = Ut/mag(Ut);

                // Other tangential unit vector
                vector tw2 = nw^tw1;

                U =
                    sqrt(physicoChemical::k.value()*T/mass)
                   *(
                        rndGen.GaussNormal<scalar>()*tw1
                      + rndGen.GaussNormal<scalar>()*tw2
                      - sqrt
                        (
                            -2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL))
                        )*nw
                    );

                scalar Un = U & nF;

                U -= 2.0*Un*nF;

                U += velocity_;

                td.switchProcessor = false;

                ++nReflections_;
            }
            else
            {
                ++nRejections_;
            }
        }
    }
    else if (d < 0) // cyclic (non-processor boundary)
    {
        if(fB != -1)
        {
            scalar pRandom = rndGen.sample01<scalar>();

            if( pRandom <= p_ ) // reflect molecule
            {
                while (mag(Ut) < SMALL)
                {
                    // If the incident velocity is parallel to the face normal, no
                    // tangential direction can be chosen.  Add a perturbation to the
                    // incoming velocity and recalculate.

                    U = vector
                    (
                        U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                        U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                        U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
                    );

                    U_dot_nw = U & nw;

                    Ut = U - U_dot_nw*nw;
                }

                // Wall tangential unit vector
                vector tw1 = Ut/mag(Ut);

                // Other tangential unit vector
                vector tw2 = nw^tw1;

                U =
                    sqrt(physicoChemical::k.value()*T/mass)
                *(
                        rndGen.GaussNormal<scalar>()*tw1
                    + rndGen.GaussNormal<scalar>()*tw2
                    - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
                    );

                scalar Un = U & nF;

                U -= 2.0*Un*nF;

                U += velocity_;

                ++nReflections_;
            }
            else
            {
                ++nRejections_;
            }
        }
    }
}


void Foam::uniGasReflectiveParticleMembranePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasReflectiveParticleMembranePatch::updateProperties
(
    const dictionary& dict
)
{
    // the main properties should be updated first
    uniGasCyclicBoundary::updateProperties(dict);
}


// ************************************************************************* //
