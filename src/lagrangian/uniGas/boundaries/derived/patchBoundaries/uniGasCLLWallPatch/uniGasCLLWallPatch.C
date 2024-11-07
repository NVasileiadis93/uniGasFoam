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

#include "uniGasCLLWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasCLLWallPatch, 0);
addToRunTimeSelectionTable(uniGasPatchBoundary, uniGasCLLWallPatch, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasCLLWallPatch::uniGasCLLWallPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    normalAccommCoeff_(propsDict_.get<scalar>("normalAccommCoeff")),
    tangentialAccommCoeff_(propsDict_.get<scalar>("tangentialAccommCoeff")),
    rotEnergyAccommCoeff_(propsDict_.get<scalar>("rotEnergyAccommCoeff")),
    temperature_(propsDict_.get<scalar>("temperature")),
    velocity_(propsDict_.get<vector>("velocity"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    Info << temperature_ << " " << velocity_.x() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasCLLWallPatch::initialConfiguration()
{
    if ((normalAccommCoeff_ < VSMALL) && (tangentialAccommCoeff_ < VSMALL))
    {
        // reduces to a specular wall, so no need to measure properties
        measurePropertiesAtWall_ = false;
    }
}


void Foam::uniGasCLLWallPatch::calculateProperties()
{}


void Foam::uniGasCLLWallPatch::controlParticle
(
    uniGasParcel& p,
    uniGasParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    label typeId = p.typeId();

    scalar& ERot = p.ERot();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

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

    const scalar T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    label rotationalDof = cloud_.constProps(typeId).rotationalDoF();

    const scalar alphaT = tangentialAccommCoeff_*(2.0 - tangentialAccommCoeff_);

    const scalar alphaN = normalAccommCoeff_;

    const scalar alphaR = rotEnergyAccommCoeff_;

    scalar mostProbableVelocity = sqrt(2.0*physicoChemical::k.value()*T/mass);

    // normalising the incident velocities

    vector normalisedTangentialVelocity = Ut/mostProbableVelocity;

    scalar normalisedNormalVelocity = U_dot_nw/mostProbableVelocity;

    // normal random number components

    scalar thetaNormal = 2.0*pi*rndGen.sample01<scalar>();

    scalar rNormal = sqrt(-alphaN*log(rndGen.sample01<scalar>()));

    // tangential random number components

    scalar thetaTangential1 = 2.0*pi*rndGen.sample01<scalar>();

    scalar rTangential1 = sqrt(-alphaT*log(rndGen.sample01<scalar>()));

    // selecting the reflected thermal velocities

    scalar normalisedIncidentTangentialVelocity1 =
        mag(normalisedTangentialVelocity);

    scalar um = sqrt(1.0 - alphaN)*normalisedNormalVelocity;

    scalar normalVelocity =
        sqrt
        (
            (rNormal*rNormal)
          + (um*um)
          + 2.0*rNormal*um*cos(thetaNormal)
        );

    scalar tangentialVelocity1 =
        sqrt(1.0 - alphaT)
       *mag(normalisedIncidentTangentialVelocity1)
      + rTangential1*cos(thetaTangential1);

    scalar tangentialVelocity2 = rTangential1*sin(thetaTangential1);

    U =
        mostProbableVelocity
       *(
            tangentialVelocity1*tw1
          + tangentialVelocity2*tw2
          - normalVelocity*nw
        );

    vector uWallNormal = (velocity_ & nw) * nw;
    vector uWallTangential1 = (velocity_ & tw1) * tw1;
    vector uWallTangential2 = (velocity_ & tw2) * tw2;
    vector UNormal = ((U & nw) * nw) + uWallNormal*alphaN;
    vector UTangential1 = (U & tw1) * tw1 + uWallTangential1*alphaT;
    vector UTangential2 = (U & tw2) * tw2 + uWallTangential2*alphaT;

    U = UNormal + UTangential1 + UTangential2;

    // selecting rotational energy, this is Lord's extension to
    // rotational degrees of freedom

    if (rotationalDof == 2)
    {
        scalar om =
            sqrt
            (
                (ERot*(1.0 - alphaR))
               /(physicoChemical::k.value()*T)
            );
        scalar rRot =
            sqrt
            (
               -alphaR*(log(max(1.0 - rndGen.sample01<scalar>(), VSMALL)))
            );
        scalar thetaRot = 2.0*pi*rndGen.sample01<scalar>();
        ERot =
            physicoChemical::k.value()*T
           *(
                (rRot*rRot)
              + (om*om)
              + (2.0*rRot*om*cos(thetaRot))
            );
    }

    if (rotationalDof == 3) // polyatomic case, see Bird's DSMC2.FOR code
    {
        scalar X = 0.0;
        scalar A = 0.0;

        do
        {
            X = 4.0*rndGen.sample01<scalar>();
            A = 2.7182818*X*X*exp(-(X*X));
        } while (A < rndGen.sample01<scalar>());

        scalar om =
            sqrt
            (
                (ERot*(1.0 - alphaR))
               /(physicoChemical::k.value()*T)
            );
        scalar rRot = sqrt(-alphaR)*X;
        scalar thetaRot = 2.0*rndGen.sample01<scalar>() - 1.0;
        ERot =
            physicoChemical::k.value()*T
           *(
                (rRot*rRot)
              + (om*om)
              + (2.0*rRot*om*cos(thetaRot))
            );
    }

    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uniGasCLLWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasCLLWallPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");
    
    velocity_ = propsDict_.get<vector>("velocity");
    
    temperature_ = propsDict_.get<scalar>("temperature");

    Info << temperature_ << " " << velocity_.x() << endl;

}

// ************************************************************************* //
