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

#include "uspLiouFangTotalPressureInlet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uspLiouFangTotalPressureInlet, 0);

addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspLiouFangTotalPressureInlet,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uspLiouFangTotalPressureInlet::uspLiouFangTotalPressureInlet
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    theta_(),
    totalPressure_(),
    totalTemperature_(),
    molecularMass_(),
    gamma_(),
    moleFractions_(),
    inletPressure_(faces_.size(), Zero),
    inletNumberDensity_(faces_.size(), Zero),
    inletTemperature_(faces_.size(), Zero),
    inletVelocity_(faces_.size(), Zero),
    previousInletPressure_(faces_.size(), Zero), 
    previousInletTemperature_(faces_.size(), Zero),       
    previousInletVelocity_(faces_.size(), Zero)
{

    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uspLiouFangTotalPressureInlet::initialConfiguration()
{}


void uspLiouFangTotalPressureInlet::calculateProperties()
{}


void uspLiouFangTotalPressureInlet::controlParcelsBeforeMove()
{
    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );
}


void uspLiouFangTotalPressureInlet::controlParcelsBeforeCollisions()
{}


void uspLiouFangTotalPressureInlet::controlParcelsAfterCollisions()
{

    step_++;

    molecularMass_ = 0.0;
    gamma_ = 0.0;
    forAll(moleFractions_, i)
    {
        const label typeId = typeIds_[i];

        molecularMass_ += cloud_.constProps(typeId).mass()*moleFractions_[i];
        gamma_ += (5.0 + cloud_.constProps(typeId).rotationalDoF())/(3.0 + cloud_.constProps(typeId).rotationalDoF())*moleFractions_[i];

    }

    scalar nParcels = 0.0;
    scalar mass = 0.0;
    scalar transT;
    scalar psi = 0.0;
    vector momentum = vector::zero;
    vector newInletVelocity= vector::zero;
    vector velSqrMean = vector::zero;
    vector velMeanSqr = vector::zero;


    const List<DynamicList<uspParcel*>>& cellOccupancy =
        cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const List<uspParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        forAll(parcelsInCell, pIC)
        {
            uspParcel* p = parcelsInCell[pIC];

            const scalar m =
                cloud_.nParticle()*cloud_.constProps(p->typeId()).mass();

            scalar RWF = cloud_.axiRWF(p->position());

            nParcels += 1.0;
            momentum += RWF*m*p->U();
            mass += RWF*m;
            velSqrMean += cmptSqr(p->U());
            velMeanSqr += p->U();
        }

        if (nParcels > 0)
        {

            velSqrMean = velSqrMean/nParcels;
            
            velMeanSqr = cmptSqr(velMeanSqr/nParcels);

            newInletVelocity = momentum/mass;

            transT =
                molecularMass_/(3.0*physicoChemical::k.value())
               *(
                    cmptSum(velSqrMean)
                  - cmptSum(velMeanSqr)
                );

            psi = molecularMass_/(gamma_*physicoChemical::k.value()*transT);

            if (step_ == 1)
            {
                previousInletPressure_ = molecularMass_*mass/(mesh_.V()[c]*physicoChemical::k.value()*transT);

                previousInletTemperature_ = transT;

                previousInletVelocity_ = newInletVelocity;                
            }

        }
        else
        {
            FatalErrorInFunction
                << "Zero parcels found at boundary cell. Decrease number of equivalent particles. "
                << abort(FatalError);
        }

        vector sF = mesh_.faceAreas()[faces_[c]];
        sF /= mag(sF);

        inletVelocity_[c] = theta_*newInletVelocity + (1.0 - theta_)*previousInletVelocity_[c];

        inletPressure_[c] = theta_*(totalPressure_-0.5*(mass/mesh_.V()[c])*pow(inletVelocity_[c] & sF, 2)) + (1.0 - theta_)*previousInletPressure_[c];

        inletTemperature_[c] = theta_*totalTemperature_/(1.0+0.5*psi*(gamma_-1.0)/gamma_*pow(inletVelocity_[c] & sF, 2)) + (1.0 - theta_)*previousInletTemperature_[c];

        if (inletPressure_[c] <= 0)
        {
            inletPressure_[c] = totalPressure_;
            inletTemperature_[c] = totalTemperature_;
        }

        inletTemperature_[c] = totalTemperature_;
        inletNumberDensity_[c] = inletPressure_[c]/(physicoChemical::k.value()*inletTemperature_[c]);

        //Info << "cellI " << c << " " << inletPressure_[c] << " " << inletTemperature_[c] << " " << inletVelocity_[c] << endl;

        std::ofstream debugFile;
        debugFile.open("debug"+std::to_string(c)+".txt", ios_base::app);
        debugFile << inletPressure_[c] << " " << inletNumberDensity_[c] << " " << mass/mesh_.V()[c] << " " << inletTemperature_[c] << " " << transT << " " << mag(inletVelocity_[c]) << nl;
        debugFile.close();

        inletTemperature_[c] = totalTemperature_;

    }

    previousInletPressure_ = inletPressure_;
    
    previousInletTemperature_ = inletTemperature_;
    
    previousInletVelocity_ = inletVelocity_;

    // Compute number of parcels to insert
    computeParcelsToInsert
    (
        inletTemperature_,
        inletVelocity_,
        inletNumberDensity_,
        moleFractions_
    );
}


void uspLiouFangTotalPressureInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void uspLiouFangTotalPressureInlet::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);

    setProperties();
}


void uspLiouFangTotalPressureInlet::setProperties()
{

    step_ = 0;

    totalPressure_ = propsDict_.get<scalar>("totalPressure");

    totalTemperature_ = propsDict_.get<scalar>("totalTemperature");

    theta_ = propsDict_.get<scalar>("theta");

    if (0.0 > theta_ || theta_ > 1.0)
    {
        FatalErrorInFunction
            << "Theta must be a value between 0 and 1 " << nl << "in: "
            << mesh_.time().system()/uspBoundaries::dictName
            << exit(FatalError);
    }

    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    // Set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }
}

} // End namespace Foam

// ************************************************************************* //
