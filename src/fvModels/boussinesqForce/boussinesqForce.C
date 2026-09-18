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

#include "boussinesqForce.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(boussinesqForce, 0);
    addToRunTimeSelectionTable(fvModel, boussinesqForce, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::boussinesqForce::boussinesqForce
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
    T_(mesh.lookupObject<volScalarField>("T")),
    thermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, alpha_.group())
        )
    ),
    beta_("beta", dimless/dimTemperature, dict.lookup<scalar>("beta")),
    T0_("T0", dimTemperature, dict.lookup<scalar>("T0")),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    if (!thermo_.isochoric())
    {
        FatalIOErrorInFunction(dict)
            << "boussinesqForce requires a constant-density equation of state "
            << "for phase " << alpha_.group() << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::boussinesqForce::addSupFields() const
{
    return wordList(1, "U");
}


void Foam::fv::boussinesqForce::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    eqn -= alpha_*thermo_.rho()*g_*beta_*(T_ - T0_);
}


void Foam::fv::boussinesqForce::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::boussinesqForce::mapMesh(const polyMeshMap&)
{}


void Foam::fv::boussinesqForce::distribute(const polyDistributionMap&)
{}


bool Foam::fv::boussinesqForce::movePoints()
{
    return true;
}


// ************************************************************************* //
