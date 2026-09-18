/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "linearExplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearExplicit, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearExplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearExplicit::linearExplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    solidificationModel(mesh, group),

    dict_(subDict(dictName_)),

    Lm_
    (
        "L",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("L")
    ),

    Tliq_
    (
        "Tliq",
        dimTemperature,
        dict_.lookup<scalar>("Tliq")
    ),

    Tsol_
    (
        "Tsol",
        dimTemperature,
        dict_.lookup<scalar>("Tsol")
    ),

    Cu_
    (
        "Cu",
        dimDensity/dimTime,
        dict_.lookupOrDefault<scalar>("Cu", 1e5)
    ),

    q_(dict_.lookupOrDefault<scalar>("q", 0.001)),

    relax_(dict_.lookupOrDefault<scalar>("relax", 1.0)),

    alphaSolid_
    (
        IOobject
        (
            IOobject::groupName("alpha.solid", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    )
{
    if (solidificationModel::debug)
    {
        const dimensionedScalar minSt
            = gMin((thermo_.Cp()*(Tliq_ - Tsol_)/Lm_)().primitiveField());
        Info<< "Minimum Stefan number St = " << minSt.value() << endl;
    }

    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::solidificationModels::linearExplicit::correct(const bool relax)
{
    sf_.storePrevIter();
    const volScalarField& sf0 = sf_.prevIter();

    // Relaxation
    if (relax)
    {
        sf_ = sf_ + relax_*thermo_.Cp()/Lm_*(Tliq_-T_ - sf_*(Tliq_-Tsol_));
        sf_ = max(min(sf_, 1.0), 0.0);
    }
    else
    {
        sf_ = max(min((Tliq_ - T_)/(Tliq_ - Tsol_), 1.0), 0.0);
    }
    alphaSolid_ = alpha_*sf_;

    if (debug)
    {
        const volScalarField sfEq =
            max(min((Tliq_ - T_)/(Tliq_ - Tsol_), 1.0), 0.0);

        Info<< "Difference to equilibrium res = "
        << gMax(mag(sf_.v() - sfEq.v())().primitiveField()) << endl;
    }

    // Residual
    return gMax(mag(sf_.v() - sf0.v())().primitiveField());
}


void Foam::solidificationModels::linearExplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    eqn -= Lm_*
        (
            fvc::ddt(alpha_, thermo_.rho(), sf_)
          + fvc::div(alphaRhoPhi_, sf_)
        );
}


void Foam::solidificationModels::linearExplicit::addSup
(
    fvMatrix<vector>& eqn
) const
{
    eqn += fvm::Sp
        (
            Cu_*sqr(alphaSolid_)/(pow3(1.0 - alphaSolid_) + q_),
            eqn.psi()
        );
}

// ************************************************************************* //
