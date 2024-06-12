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

\*----------------------------------------------------------------------------*/

#include "boundaryMeasurements.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    uniGasCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


Foam::boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const bool dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(identity(cloud_.typeIdList().size())),
    rhoNIntBF_(),
    rhoNElecBF_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    mccSpeciesBF_(),
    momentumBF_(),
    UMeanBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    qBF_(),
    fDBF_()
{
    rhoNIntBF_.setSize(typeIds_.size());
    rhoNElecBF_.setSize(typeIds_.size());
    rhoNBF_.setSize(typeIds_.size());
    rhoMBF_.setSize(typeIds_.size());
    linearKEBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    momentumBF_.setSize(typeIds_.size());
    UMeanBF_.setSize(typeIds_.size());
    rotationalEBF_.setSize(typeIds_.size());
    rotationalDofBF_.setSize(typeIds_.size());
    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    qBF_.setSize(typeIds_.size());
    fDBF_.setSize(typeIds_.size());

    forAll(rhoNBF_, i)
    {
        rhoNIntBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNElecBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        linearKEBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalDofBF_[i].setSize(mesh_.boundaryMesh().size());
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        qBF_[i].setSize(mesh_.boundaryMesh().size());
        fDBF_[i].setSize(mesh_.boundaryMesh().size());


        forAll(rhoNBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            rhoNIntBF_[i][j].setSize(patch.size(), 0.0);
            rhoNElecBF_[i][j].setSize(patch.size(), 0.0);
            rhoNBF_[i][j].setSize(patch.size(), 0.0);
            rhoMBF_[i][j].setSize(patch.size(), 0.0);
            linearKEBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            momentumBF_[i][j].setSize(patch.size(), Zero);
            UMeanBF_[i][j].setSize(patch.size(), Zero);
            rotationalEBF_[i][j].setSize(patch.size(), 0.0);
            rotationalDofBF_[i][j].setSize(patch.size(), 0.0);
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            qBF_[i][j].setSize(patch.size(), 0.0);
            fDBF_[i][j].setSize(patch.size(), Zero);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryMeasurements::clean()
{
    // Clean geometric fields

    forAll(rhoNBF_, i)
    {
        forAll(rhoNBF_[i], j)
        {
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            linearKEBF_[i][j] = 0.0;
            mccSpeciesBF_[i][j] = 0.0;
            momentumBF_[i][j] = Zero;
            UMeanBF_[i][j] = Zero;
            rotationalEBF_[i][j] = 0.0;
            rotationalDofBF_[i][j] = 0.0;
            vibrationalEBF_[i][j] = 0.0;
            electronicEBF_[i][j] = 0.0;
            qBF_[i][j] = 0.0;
            fDBF_[i][j] = Zero;

            rhoNIntBF_[i][j] = 0.0;
            rhoNElecBF_[i][j] = 0.0;
        }
    }
}


void Foam::boundaryMeasurements::updateFields(uniGasParcel& p)
{}


// ************************************************************************* //
