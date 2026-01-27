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

#include "uniGasAdsorbingWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasAdsorbingWallPatch, 0);

addToRunTimeSelectionTable
(
    uniGasPatchBoundary,
    uniGasAdsorbingWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasAdsorbingWallPatch::uniGasAdsorbingWallPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    adsorptionProbs_(),
    temperature_(propsDict_.get<scalar>("temperature")),
    velocity_(propsDict_.get<vector>("velocity"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasAdsorbingWallPatch::initialConfiguration()
{}


void Foam::uniGasAdsorbingWallPatch::calculateProperties()
{}


void Foam::uniGasAdsorbingWallPatch::controlParticle
(
    uniGasParcel& p,
    uniGasParcel::trackingData& td
)
{

    label iD(typeIds_.find(p.typeId()));

    if(iD != -1) //- particle might be adsorbed
    {
        scalar adsorbtionProbability = adsorptionProbs_[iD];
        
        if(adsorbtionProbability > cloud_.rndGen().sample01<scalar>()) //- adsorbed
        {
            //- delete the particle
             td.keepParticle = false;
        }
        else //- diffuse reflection
        {
            measurePropertiesBeforeControl(p);
            diffuseReflection(p, temperature_, velocity_);
            measurePropertiesAfterControl(p, 0.0);
        }   
    }
    else //- otherwise, it is treated as a diffuse reflection
    {
        measurePropertiesBeforeControl(p);
        diffuseReflection(p, temperature_, velocity_);
        measurePropertiesAfterControl(p, 0.0);
    }

}


void Foam::uniGasAdsorbingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasAdsorbingWallPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    uniGasPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");
    
    temperature_ = propsDict_.get<scalar>("temperature");

    velocity_ = propsDict_.get<vector>("velocity");

    setProperties();

}

void Foam::uniGasAdsorbingWallPatch::setProperties()
{

    // read in the type ids
    const List<word> molecules (propsDict_.lookup("adsorbedIds"));

    if(molecules.size() == 0)
    {
        
        FatalErrorIn("dsmcadsorbingWallPatch::setProperties()")
            << "Cannot have zero typeIds being adsorbed." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(!moleculesReduced.found(moleculeName))
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(cloud_.typeIdList().find(moleculeName));

        if(typeId == -1)
        {
            
            FatalErrorIn("dsmcadsorbingWallPatch::setProperties()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    const dictionary& adsorptionProbabilitiesDict
    (
        propsDict_.subDict("adsorptionProbabilities")
    );
    
    adsorptionProbs_.clear();

    adsorptionProbs_.setSize(typeIds_.size(), 0.0);

    forAll(adsorptionProbs_, i)
    {
        adsorptionProbs_[i] = readScalar
        (
            adsorptionProbabilitiesDict.lookup(moleculesReduced[i])
        );

    }

}

// ************************************************************************* //
