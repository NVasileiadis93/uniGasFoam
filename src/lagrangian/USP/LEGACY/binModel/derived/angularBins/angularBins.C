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

Class
    angularBins

Description

\*----------------------------------------------------------------------------*/

#include "angularBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(angularBins, 0);

addToRunTimeSelectionTable(binModel, angularBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
angularBins::angularBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    m_(propsDict_.lookup("centrePoint")),
    rVx_(propsDict_.lookup("referenceVectorX")),
    rVy_(propsDict_.lookup("referenceVectorY")),
    Rin_(propsDict_.get<scalar>("Rin")),
    Rout_(propsDict_.get<scalar>("Rout")),
    thetaStart_(propsDict_.get<scalar>("angleStart")),
    thetaEnd_(propsDict_.get<scalar>("angleEnd")),
    nBins_(propsDict_.get<scalar>("nBins")),
    area_(propsDict_.get<scalar>("thickness")),
    counterClockWise_(false)
{
    scalar deltaR = Rout_ - Rin_;

    area_ *= deltaR;

    if (deltaR <= 0.0)
    {
        FatalErrorInFunction
            << "Rout " << Rout_ << " needs to be greater than Rin: " << Rin_
            << nl << "in: "
            << "system/fieldPropertiesDict"
            << exit(FatalError);
    }

    thetaStart_ *= constant::mathematical::pi/180.0;
    thetaEnd_ *= constant::mathematical::pi/180.0;

    Info<< "thetaStart = " << thetaStart_
        << ", thetaEnd = " << thetaEnd_ << endl;

    if (thetaEnd_ <= thetaStart_)
    {
        FatalErrorInFunction
            << "angleEnd has to be larger than angleStart "
            << "- define clockwise from referenceVectorY."
            << nl << "in: "
            << "system/fieldPropertiesDict"
            << exit(FatalError);
    }

    binWidth_ = (thetaEnd_ - thetaStart_)/nBins_;

    Info << "binWidth: " << binWidth_ << endl;

    if (propsDict_.found("counterClockWise"))
    {
        counterClockWise_ = Switch(propsDict_.lookup("counterClockWise"));
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

angularBins::~angularBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// cellI is a dummy variable
label angularBins::isPointWithinBin
(
    const vector& rI,
    const label cellI
)
{
    label binNumber = -1;

    vector rIm = rI - m_;
    vector rImMod = (rVy_ & rIm)*rVy_ + (rVx_ & rIm)*rVx_;

    scalar theta = acos((rVy_ & rImMod)/mag(rImMod));
    scalar sign = rVx_ & rImMod;

    if (sign < 0.0)
    {
        theta = 2.0*constant::mathematical::pi - theta;
    }

    scalar R = mag(rImMod);

    if ((R >= Rin_) && (R <= Rout_))
    {
        if ((theta >= thetaStart_) && (theta <= thetaEnd_))
        {
            label n = label((theta - thetaStart_)/binWidth_);

            if (n >= 0)
            {
                if (n == nBins_)
                {
                    n--;
                }

                if (n < nBins_)
                {
                    binNumber = n;
                }
            }
        }
    }

    return binNumber;
}


// angles
scalarField angularBins::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}


vectorField angularBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    scalar rMean = (Rout_ + Rin_)*0.5;

    forAll(positions, i)
    {
        scalar theta = (0.5 + scalar(i))*binWidth_ + thetaStart_;
        positions[i] = m_ + rMean*cos(theta)*rVy_ + rMean*sin(theta)*rVx_;
    }

    return positions;
}


label angularBins::nBins() const
{
    return nBins_;
}


scalar angularBins::binVolume(const label n)
{
    scalar volume = area_*binWidth_*(Rout_ + Rin_)*0.5;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
