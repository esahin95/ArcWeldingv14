/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "laserDTRM.H"
#include "DTRMParticle.H"
#include "fvMatrix.H"
#include "fvcGrad.H"
#include "fvcVolumeIntegrate.H"
#include "zeroGradientFvPatchFields.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(laserDTRM, 0);
    addToRunTimeSelectionTable(fvModel, laserDTRM, dictionary);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Gather the field onto the master, concatenated in processor order
template<class Type>
void gatherAndFlatten(Foam::DynamicField<Type>& field)
{
    using namespace Foam;

    List<List<Type>> gatheredField(Pstream::nProcs());
    gatheredField[Pstream::myProcNo()] = field;
    Pstream::gatherList(gatheredField);

    field =
        ListListOps::combine<List<Type>>
        (
            gatheredField,
            accessOp<List<Type>>()
        );
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::laserDTRM::laserDTRM
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", dict.lookup<word>("phase"))
        )
    ),
    powerModelPtr_(powerModel::New(dict.subDict("powerModel"), mesh)),
    curTimeIndex_(-1),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    a_(dict.lookupOrDefault<scalar>("absorption", 1e+6)),
    reflectionModelPtr_
    (
        reflectionModel::New(dict.subDict("reflectionModel"), mesh)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    formatterPtr_(setWriter::New(dict.lookup("setFormat"), dict)),
    outputPath_
    (
        mesh.time().globalPath()/functionObjects::writeFile::outputPrefix/name()
    ),
    allPositions_(),
    allTracks_(),
    allPowers_(),
    writeIndex_(-1)
{
    mkDir(outputPath_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::laserDTRM::addSupFields() const
{
    return wordList(1, "T");
}


void Foam::fv::laserDTRM::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    Q_.storePrevIter();
    Q_ = Zero;

    const meshSearch& searchEngine = meshSearch::New(mesh());

    lagrangian::Cloud<DTRMParticle> cloud
    (
        mesh(),
        "DTRMCloud",
        IDLList<DTRMParticle>()
    );
    DTRMParticle::nTracks = 0;
    DTRMParticle::qLost = 0;

    powerModelPtr_->initialise();
    const List<vector>& positions = powerModelPtr_->positions();
    const List<scalar>& powers = powerModelPtr_->powers();
    const vector& normal = powerModelPtr_->normal();

    const scalar Qtot = sum(powers);
    DebugInfo<< "Total power initialised: Q = " << Qtot << endl;

    // Start each ray on the highest-numbered processor containing its start
    // point
    labelList cells(positions.size());
    labelList owners(positions.size());
    forAll(positions, i)
    {
        cells[i] = searchEngine.findCell(positions[i]);
        owners[i] = cells[i] != -1 ? Pstream::myProcNo() : -1;
    }
    Pstream::listCombineGather(owners, maxEqOp<label>());
    Pstream::listCombineScatter(owners);

    label nLocateBoundaryHits = 0;
    forAll(positions, i)
    {
        if (owners[i] == Pstream::myProcNo())
        {
            cloud.addParticle
            (
                new DTRMParticle
                (
                    searchEngine,
                    positions[i],
                    cells[i],
                    nLocateBoundaryHits,
                    normal,
                    powers[i]
                )
            );
        }
    }

    // Fields sampled along the rays
    const volVectorField gradAlpha(fvc::grad(alpha_));
    const volScalarField absorp(a_*(1 - alpha_));

    const interpolations::cellPoint<scalar> alphaInterp(alpha_);
    const interpolations::cellPoint<scalar> absorpInterp(absorp);
    const interpolations::cellPoint<vector> gradAlphaInterp(gradAlpha);

    allPositions_.clear();
    allTracks_.clear();
    allPowers_.clear();

    DTRMParticle::trackingData td
    (
        cloud,
        alphaInterp,
        absorpInterp,
        gradAlphaInterp,
        Q_,
        allPositions_,
        allTracks_,
        allPowers_,
        searchEngine,
        reflectionModelPtr_()
    );

    cloud.move(cloud, td);

    DebugInfo<< "Lost power fraction: "
        << returnReduce(DTRMParticle::qLost, sumOp<scalar>())/Qtot << endl;

    // Convert the absorbed power to a power density and relax in time
    Q_.primitiveFieldRef() /= mesh().V().primitiveField();
    Q_ = Q_.prevIter()*(1.0 - relax_) + Q_*relax_;

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::laserDTRM::addSup
(
    const volScalarField& rho,
    const volScalarField& T,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << nl
            << "Adding laser deposition to source: "
            << fvc::domainIntegrate(Q_) << endl;
    }

    eqn += Q_;
}


void Foam::fv::laserDTRM::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::laserDTRM::mapMesh(const polyMeshMap&)
{}


void Foam::fv::laserDTRM::distribute(const polyDistributionMap&)
{}


bool Foam::fv::laserDTRM::movePoints()
{
    return true;
}


bool Foam::fv::laserDTRM::write(const bool write) const
{
    if (Pstream::parRun())
    {
        gatherAndFlatten(allPositions_);
        gatherAndFlatten(allTracks_);
        gatherAndFlatten(allPowers_);
    }

    if (Pstream::master() && allPositions_.size())
    {
        DebugInfo<< "Writing out rays for time " << mesh().time().name()
            << " in directory " << outputPath_ << endl;

        formatterPtr_->write
        (
            outputPath_,
            IOobject::groupName("traced", Time::timeName(++writeIndex_)),
            coordSet(allTracks_, word::null, allPositions_),
            "Power",
            allPowers_
        );
    }

    return true;
}


// ************************************************************************* //
