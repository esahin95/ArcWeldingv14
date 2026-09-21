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

#include "linearSemiImplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearSemiImplicit, 0);
    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearSemiImplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidificationModels::linearSemiImplicit::linearise
(
    const scalar T,
    const scalar s,
    scalar& dsdT,
    scalar& T0
) const
{
    // Slope driving the solid fraction back to its bounds outside the
    // mushy zone
    const scalar slope = 1e10;
    const scalar tol = 1e-3;

    const scalar Tliq = Tliq_.value();
    const scalar Tsol = Tsol_.value();

    if (T > Tliq)
    {
        if (s > tol)
        {
            dsdT = slope;
            T0 = Tliq - tol;
        }
        else
        {
            dsdT = 0;
            T0 = T;
        }
    }
    else if (T < Tsol)
    {
        if (s < 1.0 - tol)
        {
            dsdT = slope;
            T0 = Tsol + tol;
        }
        else
        {
            dsdT = 0;
            T0 = T;
        }
    }
    else
    {
        dsdT = -1.0/(Tliq - Tsol);
        T0 = Tsol + (s - 1.0)/dsdT;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearSemiImplicit::linearSemiImplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    linearExplicit(mesh, group),
    T0_
    (
        IOobject
        (
            IOobject::groupName("T0", group),
            mesh.time().name(),
            mesh
        ),
        T_
    ),
    dsdT_
    (
        IOobject
        (
            IOobject::groupName("dsdT", group),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTemperature, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::solidificationModels::linearSemiImplicit::correct(const bool relax)
{
    const volScalarField sf0("sf0", sf_);

    sf_ = max(min(sf0 + dsdT_*(T_ - T0_), 1.0), 0.0);

    if (relax)
    {
        sf_ = relax_*sf_ + (1.0 - relax_)*sf0;
    }

    alphaSolid_ = alpha_*sf_;

    return gMax(mag(sf_.primitiveField() - sf0.primitiveField()));
}


void Foam::solidificationModels::linearSemiImplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    // Linearise the solid fraction in the cells ...
    forAll(dsdT_, celli)
    {
        linearise(T_[celli], sf_[celli], dsdT_[celli], T0_[celli]);
    }

    // ... and on the boundary faces
    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& Tp = T_.boundaryField()[patchi];
        const fvPatchScalarField& sfp = sf_.boundaryField()[patchi];
        fvPatchScalarField& dsdTp = dsdT_.boundaryFieldRef()[patchi];
        fvPatchScalarField& T0p = T0_.boundaryFieldRef()[patchi];

        forAll(Tp, facei)
        {
            linearise(Tp[facei], sfp[facei], dsdTp[facei], T0p[facei]);
        }
    }

    const dimensionedScalar rDeltaT = 1.0/mesh_.time().deltaT();

    eqn -=
    (
        alpha_*thermo_.rho()*Lm_*fvc::ddt(sf_)
      + rDeltaT*alpha_*thermo_.rho()*Lm_*(fvm::Sp(dsdT_, T_) - dsdT_*T0_)
    );
}


// ************************************************************************* //
