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

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable
(
    dsmcCollisionPartner,
    noTimeCounter,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noTimeCounter::noTimeCounter
(
    const dictionary& dict,
    const polyMesh& mesh,
    uniGasCloud& cloud
)
:
    dsmcCollisionPartner(dict, mesh, cloud),
    infoCounter_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noTimeCounter::initialConfiguration()
{}


void Foam::noTimeCounter::collide()
{
    if (!cloud_.dsmcCollision().active())
    {
        return;
    }

    const scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<uniGasParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModelId()[cell] == cloud_.dsmcCollModelId())
        {        

            const DynamicList<uniGasParcel*>& cellParcels(cellOccupancy[cell]);
            const scalar cellVolume = mesh.cellVolumes()[cell];

            label nC(cellParcels.size());

            if (nC > 1)
            {

                // Temporary storage for subCells
                scalar subCellLevelsProduct = 1.0; 
                vector dimWeight = vector::zero;
                const boolVector& solutionDimensions = cloud_.solutionDimensions(); 
                const vector& subCellLevels = cloud_.subCellLevels()[cell];

                forAll(solutionDimensions, dim)
                {
                    if (solutionDimensions[dim])
                    {
                        dimWeight[dim] = subCellLevelsProduct;
                        subCellLevelsProduct *= subCellLevels[dim];
                    }
                }

                vector minCellPoint = vector(GREAT, GREAT, GREAT);
                vector maxCellPoint = vector(-GREAT, -GREAT, -GREAT);
                const List<label>& cellNodes = mesh_.cellPoints()[cell];

                forAll(cellNodes, nodeI) 
                {
                    const point& cellPoint = mesh_.points()[cellNodes[nodeI]];
                    minCellPoint.x() = min(minCellPoint.x(),cellPoint.x());
                    minCellPoint.y() = min(minCellPoint.y(),cellPoint.y());
                    minCellPoint.z() = min(minCellPoint.z(),cellPoint.z());
                    maxCellPoint.x() = max(maxCellPoint.x(),cellPoint.x());
                    maxCellPoint.y() = max(maxCellPoint.y(),cellPoint.y());
                    maxCellPoint.z() = max(maxCellPoint.z(),cellPoint.z());                
                }

                const vector cellLength = maxCellPoint - minCellPoint;

                const label nSubCells = subCellLevels.x()*subCellLevels.y()*subCellLevels.z();

                List<DynamicList<label>> subCells(nSubCells);
                           
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // Assign particles to one of the virtual cartesian subCells

                // Inverse addressing specifying which subCell a parcel is in
                List<label> whichSubCell(cellParcels.size());

                vector dimPos = vector::zero;

                forAll(cellParcels, i)
                {
                    const uniGasParcel& p = *cellParcels[i];

                    vector relPos(p.position() - minCellPoint);

                    forAll(solutionDimensions, dim)
                    {
                        if (solutionDimensions[dim])
                        {
                            dimPos[dim] = label(subCellLevels[dim]*relPos[dim]/cellLength[dim]);
                        }
                    }

                    label subCell = dimPos & dimWeight;

                    subCells[subCell].append(i);

                    whichSubCell[i] = subCell;
                    
                }
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cell];

                scalar selectedPairs = 0.0;

                scalar CWF = cloud_.cellWF(cell);

                scalar RWF = 0.0;
                scalar nMols = 0.0;

                forAll(cellParcels, i)
                {
                    const uniGasParcel& p = *cellParcels[i];

                    RWF += cloud_.axiRWF(p.position());

                    nMols += 1.0;
                }

                RWF /= nMols;

                selectedPairs = + 0.5*nC*(nC - 1)*cloud_.nParticle()*CWF*RWF*sigmaTcRMax*deltaT/cellVolume;

                label nCandidates(selectedPairs);

                if (cloud_.rndGen().sample01<scalar>() < (selectedPairs-nCandidates))
                {
                    nCandidates++;
                }

                collisionCandidates += nCandidates;

                for (label c = 0; c < nCandidates; ++c)
                {
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // subCell candidate selection procedure

                    // Select the first collision candidate
                    label candidateP = rndGen_.position<label>(0, nC - 1);

                    // Declare the second collision candidate
                    label candidateQ = -1;

                    const List<label>& subCellPs =
                        subCells[whichSubCell[candidateP]];

                    label nSC = subCellPs.size();

                    if (nSC > 1)
                    {
                        // If there are two or more particle in a subCell, choose
                        // another from the same cell.  If the same candidate is
                        // chosen, choose again.

                        do
                        {
                            candidateQ =
                                subCellPs[rndGen_.position<label>(0, nSC - 1)];

                        } while (candidateP == candidateQ);
                    }
                    else
                    {
                        // Select a possible second collision candidate from the
                        // whole cell.  If the same candidate is chosen, choose
                        // again.

                        do
                        {
                            candidateQ = rndGen_.position<label>(0, nC - 1);

                        } while (candidateP == candidateQ);
                    }

                    uniGasParcel& parcelP = *cellParcels[candidateP];
                    uniGasParcel& parcelQ = *cellParcels[candidateQ];

                    label chargeP = -2;
                    label chargeQ = -2;

                    chargeP = cloud_.constProps(parcelP.typeId()).charge();
                    chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                    // do not allow electron - electron collisions

                    if (!(chargeP == -1 && chargeQ == -1))
                    {
                        scalar sigmaTcR = cloud_.dsmcCollision().sigmaTcR
                        (
                            parcelP,
                            parcelQ
                        );


                        // Update the maximum value of sigmaTcR stored, but use the
                        // initial value in the acceptance-rejection criteria
                        // because the number of collision candidates selected was
                        // based on this


                        if (sigmaTcR > cloud_.sigmaTcRMax()[cell])
                        {
                            cloud_.sigmaTcRMax()[cell] = sigmaTcR;
                        }

                        if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<scalar>())
                        {
                            // chemical reactions

                            // find which reaction model parcel p and q should use
                            label rMId =
                                cloud_.reactions().returnModelId(parcelP, parcelQ);

                            if (rMId != -1)
                            {
                                // try to react molecules
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );

                                // if reaction unsuccessful use conventional
                                // collision model
                                if (cloud_.reactions().reactions()[rMId]->relax())
                                {
                                    cloud_.dsmcCollision().collide
                                    (
                                        parcelP,
                                        parcelQ,
                                        cell
                                    );
                                }
                            }
                            else
                            {
                                // if reaction model not found, use conventional
                                // collision model
                                cloud_.dsmcCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cell
                                );
                            }

                            collisions++;
                        }
                    }
                }
            }

        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    infoCounter_++;

    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    DSMC Collisions                 = "
                << collisions << nl
                << endl;

            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No DSMC collisions" << endl;

            infoCounter_ = 0;
        }
    }
}


// ************************************************************************* //
