/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "uspMeshFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMeshFill, 0);

addToRunTimeSelectionTable
(
    uspConfiguration,
    uspMeshFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMeshFill::uspMeshFill(uspCloud& cloud, const dictionary& dict)
:
    uspConfiguration(cloud, dict)
{

    checkChapmanEnskog();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspMeshFill::checkChapmanEnskog()
{

    const word initialDistributionType(uspInitialiseDict_.lookup("initialDistributionType")); 

    if (initialDistributionType=="ChapmanEnskog")
    {

        const label n=50;   
        const scalar lBound=-5.0; 
        const scalar uBound=+5.0;
        const scalar negativeIntegralLimit=1e-9;

        Field<scalar> c(n);
        Field<scalar> wc(n);

        const scalar translationalTemperature(readScalar(uspInitialiseDict_.lookup("translationalTemperature")));

        const vector heatFlux(uspInitialiseDict_.lookup("heatFlux"));
        const tensor stress(uspInitialiseDict_.lookup("stress"));

        forAll(c, i)
            {
                c[i]=lBound+(i+0.5)*(uBound-lBound)/double(n);
                wc[i]=(uBound-lBound)/double(n);
            }

        const dictionary& numberDensitiesDict
        (
            uspInitialiseDict_.subDict("numberDensities")
        );

        List<word> molecules(numberDensitiesDict.toc());

        Field<scalar> numberDensities(molecules.size());

        scalar pressure=0.0;
        forAll(molecules, i)
        {
            numberDensities[i] = readScalar
            (
                numberDensitiesDict.lookup(molecules[i])
            );
            pressure += numberDensities[i]*physicoChemical::k.value()*translationalTemperature;
        }

        forAll(molecules, m)
            {

                const word& moleculeName = molecules[m];

                label typeId(cloud_.typeIdList().find(moleculeName));

                if (typeId == -1)
                {
                    FatalErrorIn("Foam::uspCloud<uspParcel>::initialise")
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const uspParcel::constantProperties& cP = cloud_.constProps(typeId);

                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        translationalTemperature,
                        cP.mass()
                    )
                );

                scalar negativeIntegral=0.0; 

                forAll(c, i)
                {
                    forAll(c, j)
                    {
                        forAll(c, k)
                        {

                            scalar gamma=1.0+(2.0/mostProbableSpeed*(heatFlux.x()*c[i]+heatFlux.y()*c[j]+heatFlux.z()*c[k])*(2.0/5.0*(pow(c[i],2)+pow(c[j],2)+pow(c[k],2))-1.0)
                                        -2.0*(stress.xy()*c[i]*c[j]+stress.xz()*c[i]*c[k]+stress.yz()*c[j]*c[k])
                                        -stress.xx()*(pow(c[i],2)-pow(c[k],2))-stress.yy()*(pow(c[j],2)-pow(c[k],2)))/pressure;  
             
                            

                            if (gamma<0.0)
                            {
                                scalar fM=exp(-(pow(c[i],2)+pow(c[j],2)+pow(c[k],2)))/pow(pi,1.5);
                                negativeIntegral-=fM*gamma*wc[i]*wc[j]*wc[k]; 
                            }

                        }
                    }
                }


                if (negativeIntegral>negativeIntegralLimit) 
                {
                    FatalErrorIn("Foam::uspCloud<uspParcel>::initialise")
                        << "Chapman-Enskog can not be applied too far from equilibrium" << nl
                        << abort(FatalError);        
                }

            }

    }

}

void Foam::uspMeshFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const word initialDistributionType
    (
        uspInitialiseDict_.get<word>("initialDistributionType")
    );

    const scalar translationalTemperature
    (
        uspInitialiseDict_.get<scalar>("translationalTemperature")
    );

    const scalar rotationalTemperature
    (
        uspInitialiseDict_.get<scalar>("rotationalTemperature")
    );

    const scalar vibrationalTemperature
    (
        uspInitialiseDict_.get<scalar>("vibrationalTemperature")
    );

    const scalar electronicTemperature
    (
        uspInitialiseDict_.get<scalar>("electronicTemperature")
    );

    const vector velocity(uspInitialiseDict_.get<vector>("velocity"));

    vector heatFlux=vector::zero;
    tensor stress=tensor::zero;
    if (initialDistributionType=="ChapmanEnskog")
    {

        heatFlux=uspInitialiseDict_.get<vector>("heatFlux");
        stress=uspInitialiseDict_.get<tensor>("stress");

    }

    const dictionary& numberDensitiesDict =
        uspInitialiseDict_.subDict("numberDensities");

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar>numberDensities(molecules.size());

    scalar totalNumberDensity  = 0.0;
    forAll(molecules, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
        totalNumberDensity += numberDensities[i];
        
    }

    numberDensities /= cloud_.nParticle();
    scalar pressure = totalNumberDensity*physicoChemical::k.value()*translationalTemperature;

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    //Compute cell weights
    if (cloud_.cellWeighted())
    {
        forAll(mesh_.cells(), celli)
        {

            label geometricDims = 0;
            forAll(mesh_.geometricD(), dim)
            {
                if (mesh_.geometricD()[dim] == 1)
                {
                    geometricDims++;
                }    
            }

            scalar RWF = cloud_.axiRWF(meshCC[celli]);
            cloud_.cellWeightFactor().primitiveFieldRef()[celli] =
                (totalNumberDensity*meshV[celli])/(cloud_.particlesPerSubcell()*pow(cloud_.subcellLevels()[celli],geometricDims)*cloud_.nParticle()*RWF);

        }
        cloud_.cellWeightFactor().correctBoundaryConditions();
    }
    else
    {
        cloud_.cellWeightFactor().primitiveFieldRef() = 1.0;
        cloud_.cellWeightFactor().correctBoundaryConditions();
    }

    forAll(mesh_.cells(), celli)
    {
        List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

        for (const tetIndices& cellTetIs : cellTets)
        {
            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName = molecules[i];

                label typeId(cloud_.typeIdList().find(moleculeName));

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined."
                        << abort(FatalError);
                }

                const auto& cP = cloud_.constProps(typeId);

                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar CWF = cloud_.cellWF(celli);
                scalar RWF = cloud_.axiRWF(meshCC[celli]);
                scalar particlesRequired = numberDensity*tetVolume/(CWF*RWF);

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    ++nParticlesToInsert;
                }

                // Compute breakdown parameter if particles are initialized
                //based on the Chapman-Enskog distribution
                scalar breakdownParameter=0.0;
    
                if (initialDistributionType=="ChapmanEnskog")
                {
                    scalar mostProbableSpeed(
                        cloud_.maxwellianMostProbableSpeed
                        (
                            translationalTemperature,
                            cP.mass()
                        )
                    );

                    scalar maxHeatFlux=-1.0;
                    forAll(heatFlux,i) 
                    {
                        if (maxHeatFlux < fabs(heatFlux[i])) 
                        {
                            maxHeatFlux = fabs(heatFlux[i]);
                        }
                    }
                    maxHeatFlux = 2.0*maxHeatFlux/(pressure*mostProbableSpeed);

                    scalar maxStress=-1.0;
                    forAll(stress,i) 
                    {
                        if (maxStress < fabs(stress[i])) 
                        {
                            maxStress = fabs(stress[i]);
                        }
                    }
                    maxStress = maxStress/pressure;

                    breakdownParameter = max(maxHeatFlux,maxStress);

                }

                for (label pI = 0; pI < nParticlesToInsert; ++pI)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U=vector::zero;
                    if (initialDistributionType=="Maxwellian")
                    {

                        U = cloud_.equipartitionLinearVelocity
                        (
                            translationalTemperature,
                            cP.mass()
                        );
                        
                    }
                    else if (initialDistributionType=="ChapmanEnskog")
                    {
                        U = cloud_.chapmanEnskogLinearVelocity
                        (
                            breakdownParameter,
                            pressure,
                            translationalTemperature,
                            cP.mass(),
                            heatFlux,
                            stress
                        );

                    }
                    else
                    {
                        FatalErrorIn("Foam::uspCloud<uspParcel>::initialise")
                            << "initialDistributionType " << initialDistributionType << " not defined." << nl
                            << abort(FatalError);
                    }


                    scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        rotationalTemperature,
                        cP.rotationalDoF()
                    );

                    labelList vibLevel =
                        cloud_.equipartitionVibrationalEnergyLevel
                        (
                            vibrationalTemperature,
                            cP.vibrationalDoF(),
                            typeId
                        );

                    label ELevel = cloud_.equipartitionElectronicLevel
                    (
                        electronicTemperature,
                        cP.degeneracyList(),
                        cP.electronicEnergyList(),
                        typeId
                    );

                    U += velocity;

                    label newParcel = 0;

                    scalar CWF = cloud_.cellWF(celli);
                    scalar RWF = cloud_.axiRWF(meshCC[celli]);

                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        CWF,
                        RWF,
                        ERot,
                        ELevel,
                        celli,
                        typeId,
                        newParcel,
                        vibLevel
                    );
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222 - 223)

    label mostAbundantType(findMax(numberDensities));

    const auto& cP = cloud_.constProps(mostAbundantType);

    cloud_.sigmaTcRMax().primitiveFieldRef() =
        cP.sigmaT()
       *cloud_.maxwellianMostProbableSpeed
        (
            translationalTemperature,
            cP.mass()
        );

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


// ************************************************************************* //
