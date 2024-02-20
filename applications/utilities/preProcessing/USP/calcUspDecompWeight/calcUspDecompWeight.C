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

Application
    calcUspDecompWeight

Description
    Calculate the load balancing weight from the initial number density and cell volumes

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOdictionary uspProperties
    (
        IOobject
        (
            "uspProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read in the type ids
    List<word> molecules = uspProperties.lookup("typeIdList");

    //- List of inlet densities (one entry for each species)
    List<autoPtr<volScalarField>> initialNumberDensityPtr;

    initialNumberDensityPtr.setSize(molecules.size());

    forAll(molecules, m)
    {

        const word& moleculeName = molecules[m];

        initialNumberDensityPtr[m].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "numberDensity_"+moleculeName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

    }

    // Read number of equivalent particles
    scalar nParticle = uspProperties.get<scalar>("nEquivalentParticles");

    Switch cellWeighted(uspProperties.get<Switch>("cellWeightedSimulation"));
    scalar particlesPerCell = 0.0;
    if (cellWeighted)
    {
        particlesPerCell = uspProperties.get<scalar>("particlesPerCell");
    }

    // Read if simulation is axisymmetric
    Switch axisymmetric(uspProperties.get<Switch>("axisymmetricSimulation"));
    scalar radialExtent = 0.0;
    scalar maxRWF = 0.0;
    if (axisymmetric)
    {
        radialExtent = uspProperties.get<scalar>("radialExtentOfDomain");
        maxRWF = uspProperties.get<scalar>("maxRadialWeightingFactor");
    }

    // Perform calculation for 0 directory
    instantList timeDirs = timeSelector::select0(runTime, args);

    runTime.setTime(timeDirs[0], 0);

    Info << "Calculating usp decomposition weight field for time " << runTime.timeName() << nl << endl;

    volScalarField uspDecompositionWeight
    (
        IOobject
        (
            "uspDecompositionWeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    // Calculate usp particles to be inserted in each cell
    forAll(mesh.cells(), i)
    {

        if (cellWeighted)
        {
            uspDecompositionWeight[i] = particlesPerCell;
        }
        else
        {

            scalar totalNumberDensity = 0.0;
            forAll(molecules, m)
            {
                const volScalarField& initialNumberDensity = initialNumberDensityPtr[m];
                totalNumberDensity += initialNumberDensity[i];
            }
            totalNumberDensity /= nParticle;

            scalar RWF = 1.0;
            if (axisymmetric)
            {
                RWF = 1.0 + (maxRWF-1.0)*mesh.C()[i].y()/radialExtent;
            } 
            uspDecompositionWeight[i] = totalNumberDensity*mesh.V()[i]/RWF;

        }

    }

    // Output decomposition weight
    uspDecompositionWeight.write();

    Info << "End" << endl;

    return(0);

}


// ************************************************************************* //
