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

#include "uspPatchBoundary.H"
#include "graph.H"
#include "uspCloud.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspPatchBoundary, 0);
defineRunTimeSelectionTable(uspPatchBoundary, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspPatchBoundary::uspPatchBoundary
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspBoundaryBase(mesh, cloud, dict, "patch"),
    faces_(),
    patchSurfaceArea_(Zero),
    totalPatchSurfaceArea_(Zero),
    cells_(),
    density_(Zero),
    velocity_(Zero),
    temperature_(Zero),
    measurePropertiesAtWall_(false),
    preIE_(Zero),
    preIMom_(Zero)
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    // Initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    // Loop through all faces and set the boundary cells
    // - no conflict with parallelisation because the faces are unique

    if (isA<cyclicPolyPatch>(patch))
    {
        FatalErrorInFunction
            << "Patch: " << patchName_ << " is a cyclic boundary. "
            << "It should be a patch."
            << nl << "in: "
            << mesh_.time().system()/uspBoundaries::dictName
            << exit(FatalError);
    }

    for (label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    totalPatchSurfaceArea_ = patchSurfaceArea_;

    if (Pstream::parRun())
    {
        // Total area on all processors
        reduce(totalPatchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uspPatchBoundary> Foam::uspPatchBoundary::New
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
{
    word modelType(dict.get<word>("boundaryModel"));

    Info<< "Selecting uspPatchBoundaryModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uspPatchBoundary",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspPatchBoundary>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspPatchBoundary::measurePropertiesBeforeControl
(
    const uspParcel& p
)
{
    if (measurePropertiesAtWall_)
    {
        const label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());

        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

        const label typei = p.typeId();
        const auto& constProps = cloud_.constProps(typei);

        const scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        scalar U_dot_nw = p.U() & nw;

        const vector Ut = p.U() - U_dot_nw*nw;

        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);

        boundaryMeasurements& bm = cloud_.boundaryFluxMeasurements();
        bm.rhoNBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;

        if (constProps.rotationalDoF() > 0)
        {
            bm.rhoNIntBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;
        }

        if (constProps.nElectronicLevels() > 1)
        {
            bm.rhoNElecBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;
        }

        bm.rhoMBF()[typei][wppIndex][wppLocalFace] += m*invMagUnfA;
        bm.linearKEBF()[typei][wppIndex][wppLocalFace] +=
            0.5*m*(p.U() & p.U())*invMagUnfA;
        bm.mccSpeciesBF()[typei][wppIndex][wppLocalFace] +=
            m*(p.U() & p.U())*invMagUnfA;
        bm.momentumBF()[typei][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
        bm.rotationalEBF()[typei][wppIndex][wppLocalFace] +=
            p.ERot()*invMagUnfA;
        bm.rotationalDofBF()[typei][wppIndex][wppLocalFace] +=
            constProps.rotationalDoF()*invMagUnfA;

        if (constProps.vibrationalDoF() > 0)
        {
            forAll(p.vibLevel(), i)
            {
                bm.vibrationalEBF()[typei][wppIndex][wppLocalFace] +=
                    p.vibLevel()[i]
                   *constProps.thetaV()[i]
                   *physicoChemical::k.value()*invMagUnfA;
            }
        }

        bm.electronicEBF()[typei][wppIndex][wppLocalFace] +=
            constProps.electronicEnergyList()[p.ELevel()]*invMagUnfA;

        // Pre-interaction energy
        preIE_ =
            0.5*m*(p.U() & p.U())
          + p.ERot()
          + constProps.electronicEnergyList()[p.ELevel()];

        if (constProps.vibrationalDoF() > 0)
        {
            forAll(p.vibLevel(), i)
            {
                preIE_ +=
                    p.vibLevel()[i]
                   *constProps.thetaV()[i]
                   *physicoChemical::k.value();
            }
        }

        // pre-interaction momentum
        preIMom_ = m*p.U();
    }
}


void Foam::uspPatchBoundary::measurePropertiesAfterControl
(
    const uspParcel& p,
    const scalar heatOfReaction
)
{
    if (measurePropertiesAtWall_)
    {
        const label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());
        const label wppCell = p.cell();
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

        const scalar deltaT = mesh_.time().deltaTValue();

        const label typei = p.typeId();
        const uspParcel::constantProperties& constProps
        (
            cloud_.constProps(typei)
        );

        const scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        scalar U_dot_nw = p.U() & nw;

        vector Ut = p.U() - U_dot_nw*nw;

        scalar invMagUnfA = 1.0/max(mag(U_dot_nw)*fA, VSMALL);

        boundaryMeasurements& bm = cloud_.boundaryFluxMeasurements();

        bm.rhoNBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;

        if (constProps.rotationalDoF() > 0)
        {
            bm.rhoNIntBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;
        }
        if (constProps.nElectronicLevels() > 1)
        {
            bm.rhoNElecBF()[typei][wppIndex][wppLocalFace] += invMagUnfA;
        }

        bm.rhoMBF()[typei][wppIndex][wppLocalFace] += m*invMagUnfA;
        bm.linearKEBF()[typei][wppIndex][wppLocalFace] +=
            0.5*m*(p.U() & p.U())*invMagUnfA;
        bm.mccSpeciesBF()[typei][wppIndex][wppLocalFace] +=
            m*(p.U() & p.U())*invMagUnfA;
        bm.momentumBF()[typei][wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
        bm.rotationalEBF()[typei][wppIndex][wppLocalFace] +=
            p.ERot()*invMagUnfA;
        bm.rotationalDofBF()[typei][wppIndex][wppLocalFace] +=
            constProps.rotationalDoF()*invMagUnfA;
        forAll(p.vibLevel(), i)
        {
            bm.vibrationalEBF()[typei][wppIndex][wppLocalFace] +=
                p.vibLevel()[i]
               *constProps.thetaV()[i]
               *physicoChemical::k.value()
               *invMagUnfA;
        }
        bm.electronicEBF()[typei][wppIndex][wppLocalFace] +=
            constProps.electronicEnergyList()[p.ELevel()]*invMagUnfA;

        // post - interaction energy
        scalar postIE =
            0.5*m*(p.U() & p.U())
          + p.ERot()
          + constProps.electronicEnergyList()[p.ELevel()];

        if (constProps.vibrationalDoF() > 0)
        {
            forAll(p.vibLevel(), i)
            {
                postIE +=
                    p.vibLevel()[i]
                   *constProps.thetaV()[i]
                   *physicoChemical::k.value();
            }
        }
        // post - interaction momentum
        const vector postIMom = m*p.U();

        scalar CWF = cloud_.cellWF(wppCell);
        scalar RWF = cloud_.axiRWF(p.position());
        scalar nParticle = cloud_.nParticle()*CWF*RWF;

        scalar deltaQ =
            nParticle
           *(preIE_ - postIE + (heatOfReaction*physicoChemical::k.value()))
           /(deltaT*fA);
        vector deltaFD = nParticle*(preIMom_ - postIMom)/(deltaT*fA);

        bm.qBF()[typei][wppIndex][wppLocalFace] += deltaQ;
        bm.fDBF()[typei][wppIndex][wppLocalFace] += deltaFD;
    }
}


void Foam::uspPatchBoundary::diffuseReflection
(
    uspParcel& p,
    const scalar wallTemp,
    const vector& wallVelocity
)
{
    vector& U = p.U();

    scalar& ERot = p.ERot();

    labelList& vibLevel = p.vibLevel();

    label& ELevel = p.ELevel();

    label typeId = p.typeId();

    vector nw(p.normal());
    nw.normalise();

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
    vector tw1(Ut/mag(Ut));

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    const scalar T = wallTemp;

    scalar mass = cloud_.constProps(typeId).mass();

    label rotationalDof = cloud_.constProps(typeId).rotationalDoF();

    label vibrationalDof = cloud_.constProps(typeId).vibrationalDoF();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal<scalar>()*tw1
          + rndGen.GaussNormal<scalar>()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
        );

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);

    vibLevel =
        cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);

    ELevel =
        cloud_.equipartitionElectronicLevel
        (
            T,
            cloud_.constProps(typeId).degeneracyList(),
            cloud_.constProps(typeId).electronicEnergyList(),
            typeId
        );

    U += wallVelocity;
}

void Foam::uspPatchBoundary::specularReflection(uspParcel& p)
{
    vector& U = p.U();

    vector nw = p.normal();
    nw.normalise();

    const scalar U_dot_nw = U & nw;

    if (U_dot_nw > 0.0)
    {
        U -= 2.0*U_dot_nw*nw;
    }
}


const Foam::labelList& Foam::uspPatchBoundary::controlPatch() const
{
    return faces_;
}


const Foam::labelList& Foam::uspPatchBoundary::controlZone() const
{
    return cells_;
}


Foam::scalar Foam::uspPatchBoundary::density() const
{
    return density_;
}


Foam::scalar& Foam::uspPatchBoundary::density()
{
    return density_;
}


const Foam::vector& Foam::uspPatchBoundary::velocity() const
{
    return velocity_;
}


Foam::vector& Foam::uspPatchBoundary::velocity()
{
    return velocity_;
}


Foam::scalar Foam::uspPatchBoundary::temperature() const
{
    return temperature_;
}


Foam::scalar& Foam::uspPatchBoundary::temperature()
{
    return temperature_;
}


// ************************************************************************* //
