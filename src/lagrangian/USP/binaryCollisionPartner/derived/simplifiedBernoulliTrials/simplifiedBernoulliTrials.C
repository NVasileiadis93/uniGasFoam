/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "simplifiedBernoulliTrials.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simplifiedBernoulliTrials, 0);

addToRunTimeSelectionTable(binaryCollisionPartner, simplifiedBernoulliTrials, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
simplifiedBernoulliTrials::simplifiedBernoulliTrials
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    binaryCollisionPartner(mesh, cloud, dict)
{ }



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simplifiedBernoulliTrials::~simplifiedBernoulliTrials()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simplifiedBernoulliTrials::initialConfiguration()
{

    Info << "Calculating simplifiedBernoulliTrials subcell volumes" << endl;

    // Monte Carlo integration to find subcell volumes
    const label integrationPoints = 10000;

    const vector solutionDims = cloud_.mesh().solutionD();

    vector dimensionMultiplier = vector::zero;
    label nSolutionDims = 0;
    forAll(solutionDims, i)
    {
        if (solutionDims[i] == 1)
        {
            dimensionMultiplier[i] = pow(2,nSolutionDims);
            nSolutionDims++;
        }
    }
    
    const label nSubcells = pow(2,nSolutionDims);
    subcellVolume = List<List<scalar>>(cloud_.mesh().nCells(),List<scalar>(nSubcells,scalar(0)));   

    forAll(mesh_.cells(), cellI)
    {

        label realIntegationPoints = 0;

        point cC = cloud_.mesh().cellCentres()[cellI];

        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices(mesh_, cellI);

        // Create random points in all cell tets
        for (const tetIndices& cellTetIs : cellTets)
        {
            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            label tetIntegrationPoints = integrationPoints*(tetVolume/cloud_.mesh().cellVolumes()[cellI]);
            realIntegationPoints += tetIntegrationPoints;

            for (label pI = 0; pI < tetIntegrationPoints; ++pI)
            {
                
                point p = tet.randomPoint(rndGen_);

                vector relPos = p - cC;

                // Find in which subcells the point is in
                label subCell = pos(relPos.x())*dimensionMultiplier.x() 
                              + pos(relPos.y())*dimensionMultiplier.y()  
                              + pos(relPos.z())*dimensionMultiplier.z();

                subcellVolume[cellI][subCell] += 1.0;

            }

        }

        // Calculate subcell volume based on subcell point and cell volume
        forAll(subcellVolume[cellI], subcellI)
        {
            subcellVolume[cellI][subcellI] *= cloud_.mesh().cellVolumes()[cellI]/scalar(realIntegationPoints);
        } 

    }

}

void simplifiedBernoulliTrials::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    const vector solutionDims = cloud_.mesh().solutionD();

    vector dimensionMultiplier = vector::zero;
    label nSolutionDims = 0;
    forAll(solutionDims, i)
    {
        if (solutionDims[i] == 1)
        {
            dimensionMultiplier[i] = pow(2,nSolutionDims);
            nSolutionDims++;
        }
    }
    const label nSubcells = pow(2,nSolutionDims);

    List<DynamicList<label>> subCells(nSubcells);

    scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisions = 0;
	
	List<DynamicList<uspParcel*>> cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {
        DynamicList<uspParcel*>& cellParcels(cellOccupancy[cellI]);

        rndGen_.shuffle(cellParcels);

        /*if (cellI == cellOccupancy.size()/2)
        {
            Info << "Cell: " << cellI << endl;
            Info << "Cell parcels:" << endl;
            forAll(cellParcelsShuffled, p)
            {
                Info << cellParcelsShuffled[p] << endl;
            }
        }

        rndGen_.shuffle(cellParcelsShuffled);

        if (cellI == cellOccupancy.size()/2)
        {
            Info << "Cell parcels shuffled:" << endl;
            forAll(cellParcelsShuffled, p)
            {
                Info << cellParcelsShuffled[p] << endl;
            }
        }*/

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            point cC = mesh.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const uspParcel& p = *cellParcels[i];

                vector relPos = p.position() - cC;

                label subCell = pos(relPos.x())*dimensionMultiplier.x() 
                              + pos(relPos.y())*dimensionMultiplier.y()  
                              + pos(relPos.z())*dimensionMultiplier.z();

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            scalar nParticle = cloud_.nParticle();
            
            scalar CWF = cloud_.cellWF(cellI);

            scalar RWF = cloud_.axiRWF(cC);
            
            nParticle *= CWF*RWF;
            
            label k = -1;
            label candidateP = -1;
            label candidateQ = -1;

            /*if (cellI == 150)
            {
            Info << "-------------------------------------" << endl;
            Info << "Cell i=" << cellI << endl;
            Info << "Number of cell particles: " << nC << endl;
            }*/

            // loop over sub cells
            forAll(subCells, i)
            {

                label nCS = subCells[i].size();
                /*if (cellI == 150)
                {
                    Info << "Subcell j=" << i << endl;
                    Info << "Number of subcell particles: " << nCS << endl;
                }*/

                if(nCS > 1)
                {
                    for(label p = 1 ; p <= nCS-1 ; p++)
                    {                
                        // Select the first collision candidate
                        candidateP = p-1;
                    
                        k = nCS - p;
                        candidateQ = candidateP + rndGen_.position<label>(1, k);

                        uspParcel& parcelP = *cellParcels[subCells[i][candidateP]];
                        uspParcel& parcelQ = *cellParcels[subCells[i][candidateQ]];  

                        scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                        (
                            parcelP,
                            parcelQ
                        );
                    
                        scalar probability = k*(nParticle*deltaT)/(cloud_.mesh().cellVolumes()[cellI]/scalar(nSubcells))*sigmaTcR;

                        /*if (cellI == 150)
                        {
                            Info << "Candidates " << candidateP << " " << candidateQ << " " << probability << endl;
                        }*/

                        if (probability > rndGen_.sample01<scalar>())
                        {
                            // chemical reactions

                            // find which reaction model parcel p and q should use
                            label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

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

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    infoCounter_++;

    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisions > 0)
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

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
