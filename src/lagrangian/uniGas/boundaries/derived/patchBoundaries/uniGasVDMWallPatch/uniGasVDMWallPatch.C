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

#include "uniGasVDMWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasVDMWallPatch, 0);

addToRunTimeSelectionTable
(
    uniGasPatchBoundary,
    uniGasVDMWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasVDMWallPatch::uniGasVDMWallPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    accomCoeff_(propsDict_.get<scalar>("accomCoeff")),
    temperature_(propsDict_.get<scalar>("temperature")),
    velocity_(propsDict_.get<vector>("velocity")),
    roughnessDegree_(propsDict_.get<scalar>("roughnessDegree")),
    activationEnergy_(propsDict_.get<scalar>("activationEnergy")),
    impactCoeff_(propsDict_.get<scalar>("impactCoeff"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasVDMWallPatch::initialConfiguration()
{}


void Foam::uniGasVDMWallPatch::calculateProperties()
{}


void Foam::uniGasVDMWallPatch::controlParticle
(
    uniGasParcel& p,
    uniGasParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();
    scalar& ERot = p.ERot();
    labelList& vibLevel = p.vibLevel();
    label& ELevel = p.ELevel();
    label typeId = p.typeId();
    scalar mass = cloud_.constProps(typeId).mass();

    vector nw(p.normal());
    nw.normalise();
    scalar U_dot_nw = U & nw;
    vector Ut = U - U_dot_nw*nw;

    while (mag(Ut) < SMALL)
    {
       	// If the incident velocity is parallel to the face normal, no
       	// tangential direction can be chosen.  Add a perturbation to the
       	// incoming velocity and recalculate.

	    U = vector
            (
     	       U.x()*(0.8 + 0.2*cloud_.rndGen().sample01<scalar>()),
               U.y()*(0.8 + 0.2*cloud_.rndGen().sample01<scalar>()),
               U.z()*(0.8 + 0.2*cloud_.rndGen().sample01<scalar>())
            );

	    U_dot_nw = U & nw;
            Ut = U - U_dot_nw*nw;
     }


    scalar phi = accomCoeff_*exp((activationEnergy_ - impactCoeff_*(mass/2)*sqr(mag(U)))/(physicoChemical::k.value()*temperature_));

    scalar R1 = cloud_.rndGen().sample01<scalar>();
    scalar R2 = cloud_.rndGen().sample01<scalar>();


    if (phi >= R1)
    {
	    // Diffuse reflection

	    scalar beta = sqrt(mass/(2*physicoChemical::k.value()*temperature_));
	    scalar Fr = cloud_.rndGen().sample01<scalar>();
	    scalar Fn = cloud_.rndGen().sample01<scalar>();
	    scalar Fo = cloud_.rndGen().sample01<scalar>();
	    U = vector
	    (
	     	(1/sqrt(1+impactCoeff_))*(1/beta)*sqrt(-log(1-Fr))*cos(2*constant::mathematical::pi*Fo),
		(1/sqrt(1+impactCoeff_))*(1/beta)*sqrt(-log(1-Fn)),
		(1/sqrt(1+impactCoeff_))*(1/beta)*sqrt(-log(1-Fr))*sin(2*constant::mathematical::pi*Fo)
	    );

	    ERot = cloud_.equipartitionRotationalEnergy(temperature_, cloud_.constProps(typeId).rotationalDoF());

    	    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(temperature_, cloud_.constProps(typeId).vibrationalDoF(), typeId);

    	    ELevel =
        	cloud_.equipartitionElectronicLevel
        	(
            		temperature_,
            		cloud_.constProps(typeId).degeneracyList(),
            		cloud_.constProps(typeId).electronicEnergyList(),
            		typeId
        	);

	    U += velocity_;
    }
    else if (accomCoeff_ >= R2)
    {
        // Specular reflection
	
        specularReflection(p);
    }
    else
    {
	    // Isotropic reflection

	    scalar beta = sqrt(mass/(2*physicoChemical::k.value()*temperature_));
	    scalar Ft = cloud_.rndGen().sample01<scalar>();
	    scalar Fp = cloud_.rndGen().sample01<scalar>();
	    U = vector
	    (
	     	mag(U)*sqrt(Fp)*cos(2*constant::mathematical::pi*Ft),
		mag(U)*sqrt(1-Fp),
		mag(U)*sqrt(Fp)*sin(2*constant::mathematical::pi*Ft)
	    );

	    ERot = cloud_.equipartitionRotationalEnergy(temperature_, cloud_.constProps(typeId).rotationalDoF());

    	    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(temperature_, cloud_.constProps(typeId).vibrationalDoF(), typeId);

    	    ELevel =
        	cloud_.equipartitionElectronicLevel
        	(
            		temperature_,
            		cloud_.constProps(typeId).degeneracyList(),
            		cloud_.constProps(typeId).electronicEnergyList(),
            		typeId
        	);

	    U += velocity_;
    }


    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uniGasVDMWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasVDMWallPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uniGasPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    temperature_ = propsDict_.get<scalar>("temperature");

    velocity_ = propsDict_.get<vector>("velocity");

}

// ************************************************************************* //
