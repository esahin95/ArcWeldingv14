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

#include "linearImplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcGrad.H"

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearImplicit, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearImplicit,
        dictionary
    );
}
}

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearImplicit::linearImplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    linearExplicit(mesh, group),

    rhoCpLatent_
    (
        IOobject
        (
            IOobject::groupName("rhoCpLatent", group),
            mesh.time().name(),
            mesh
        ),
        alpha_*thermo_.rho()*thermo_.Cp()
    ),

    rhoPhiCpLatent_
    (
        IOobject
        (
            IOobject::groupName("rhoPhiCpLatent", group),
            mesh.time().name(),
            mesh
        ),
        alphaRhoPhi_*fvc::interpolate(thermo_.Cp())
    )
{
    rhoCpLatent_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar
Foam::solidificationModels::linearImplicit::correct(const bool relax)
{
    sf_.storePrevIter();
    const volScalarField& sf0 = sf_.prevIter();

    sf_ = max(min((Tliq_ - T_)/(Tliq_ - Tsol_), 1.0), 0.0);
    if (relax)
    {
        sf_ = relax_*sf_ + (1.0 - relax_)*sf0;
    }

    alphaSolid_ = alpha_*sf_;

    // Residual
    return gMax(mag(sf_.v() - sf0.v())().primitiveField());
}


void Foam::solidificationModels::linearImplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    const scalar TOL = 1e-6;
    rhoCpLatent_ =
        Lm_/(Tliq_-Tsol_)*alpha_*thermo_.rho()*pos(sf_-TOL)*pos(1.0-TOL-sf_);
    eqn += fvm::ddt(rhoCpLatent_, T_) - fvm::Sp(fvc::ddt(rhoCpLatent_), T_);
}


// ************************************************************************* //
