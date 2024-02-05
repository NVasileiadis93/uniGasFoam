/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
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
    binaryCollisionPartner,
    noTimeCounter,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    binaryCollisionPartner(mesh, cloud, dict),
    infoCounter_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noTimeCounter::initialConfiguration()
{}


void Foam::noTimeCounter::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    const scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<uspParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {

        if (cloud_.cellCollModel(cellI) == cloud_.binCollModel())
        {        

            // Temporary storage for subCells
            label nSolutionDims = 0;
            vector solutionDims; 
            vector dimWeight = vector::zero;
            const label nSubcellLevels = cloud_.subcellLevels()[cellI];
            
            if (!cloud_.axisymmetric())
            {
                forAll(solutionDims, dim)
                {
                    solutionDims[dim] = cloud_.mesh().solutionD()[dim];
                    if (solutionDims[dim] == 1)
                    {
                        dimWeight[dim] = pow(nSubcellLevels,nSolutionDims);
                        nSolutionDims++;
                    }
                }
            }
            else
            {
                solutionDims.x() = 1;
                solutionDims.y() = 1;
                solutionDims.z() = -1;
                dimWeight.x() = 1;
                dimWeight.y() = nSubcellLevels;
                dimWeight.z() = 0;
                nSolutionDims = 2;
            }    

            vector minCellPoint = vector(GREAT, GREAT, GREAT);
            vector maxCellPoint = vector(-GREAT, -GREAT, -GREAT);
            const List<label>& cellNodes = mesh_.cellPoints()[cellI];

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
            
            const label nSubcells = pow(nSubcellLevels,nSolutionDims);

            List<DynamicList<label>> subCells(nSubcells);

            const vector cellLength = maxCellPoint - minCellPoint;

            //std::cout << "-------------------------------------------------------"  << std::endl;
            //std::cout << nSolutionDims << " " << nSubcellLevels << " " << nSubcells << std::endl;
            //std::cout << solutionDims.x() << " " << solutionDims.y() << " " << solutionDims.z() << std::endl;
            //std::cout << "-------------------------------------------------------"  << std::endl;
            //std::cout << cellLength.x() << " " << cellLength.y() << " " << cellLength.z() << std::endl;
            //std::cout << "-------------------------------------------------------"  << std::endl;
            //std::cin.get();

            const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cellI]);
            const scalar cellVolume = mesh.cellVolumes()[cellI];

            label nC(cellParcels.size());

            if (nC > 1)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // Assign particles to one of the Cartesian subCells

                // Clear temporary lists
                for (auto& sc : subCells)
                {
                    sc.clear();
                }

                // Inverse addressing specifying which subCell a parcel is in
                List<label> whichSubCell(cellParcels.size());

                vector dimPos = vector::zero;

                forAll(cellParcels, i)
                {
                    const uspParcel& p = *cellParcels[i];

                    vector relPos(p.position() - minCellPoint);

                    forAll(solutionDims, dim)
                    {
                        if (solutionDims[dim] == 1)
                        {
                            for (label l = 0; l < nSubcellLevels; ++l)
                            {
                                
                                if (relPos[dim] <= (l+1)/scalar(nSubcellLevels)*cellLength[dim])
                                {
                                    dimPos[dim] = l;
                                    break;
                                }
                            }
                        }
                    }

                    label subCell = dimPos & dimWeight;

                    subCells[subCell].append(i);

                    whichSubCell[i] = subCell;
                }

                //forAll(subCells, i)
                //{
                //    std::cout << "SC: " << i << " " << subCells[i].size() << std::endl;
                //}
                //std::cin.get();

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

                scalar selectedPairs = 0.0;

                scalar CWF = cloud_.cellWF(cellI);

                scalar RWF = 0.0;
                scalar nMols = 0.0;

                forAll(cellParcels, i)
                {
                    const uspParcel& p = *cellParcels[i];

                    RWF += cloud_.axiRWF(p.position());

                    nMols += 1.0;
                }

                RWF /= nMols;

                selectedPairs =
                    cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticle()*CWF*RWF*sigmaTcRMax*deltaT
                    /cellVolume;

                label nCandidates(selectedPairs);

                cloud_.collisionSelectionRemainder()[cellI] =
                    selectedPairs - nCandidates;

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

                    uspParcel& parcelP = *cellParcels[candidateP];
                    uspParcel& parcelQ = *cellParcels[candidateQ];

                    label chargeP = -2;
                    label chargeQ = -2;

                    chargeP = cloud_.constProps(parcelP.typeId()).charge();
                    chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                    // do not allow electron - electron collisions

                    if (!(chargeP == -1 && chargeQ == -1))
                    {
                        scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                        (
                            parcelP,
                            parcelQ
                        );


                        // Update the maximum value of sigmaTcR stored, but use the
                        // initial value in the acceptance-rejection criteria
                        // because the number of collision candidates selected was
                        // based on this


                        if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                        {
                            cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
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
                                    cloud_.binaryCollision().collide
                                    (
                                        parcelP,
                                        parcelQ,
                                        cellI
                                    );
                                }
                            }
                            else
                            {
                                // if reaction model not found, use conventional
                                // collision model
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
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
            Info<< "    Collisions                      = "
                << collisions << nl
                << endl;

            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;

            infoCounter_ = 0;
        }
    }
}


// ************************************************************************* //
