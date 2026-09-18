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
#include "fvmDiv.H"
#include "fvmSup.H"

#include "mathematicalConstants.H"

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

using namespace Foam::constant;


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
        dimensionedScalar(dimless/dimTemperature, 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar
Foam::solidificationModels::linearSemiImplicit::correct(const bool relax)
{
    const volScalarField sf0("sf0", sf_);

    sf_ = max(min(sf0 + dsdT_*(T_-T0_), 1.0), 0.0);
    if (relax)
    {
        sf_ = relax_*sf_ + (1.0-relax_)*sf0;
    }

    alphaSolid_ = alpha_*sf_;

    // Residual
    return gMax(mag(sf_.v().primitiveField() - sf0.v().primitiveField()));
}


void Foam::solidificationModels::linearSemiImplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    // Update slope
    dsdT_ = 1.0/(Tliq_ - Tsol_);
    const scalar slope = 1e10;
    const scalar tol = 1e-3;

    forAll(dsdT_, cellI)
    {
        const scalar T = T_[cellI];
        const scalar s = sf_[cellI];
        if (T > Tliq_.value())
        {
            if (s > tol)
            {
                dsdT_[cellI] = slope;
                T0_[cellI] = Tliq_.value() - tol;
            }
            else
            {
                dsdT_[cellI] = 0.0;
                T0_[cellI] = T;
            }
        }
        else if (T < Tsol_.value())
        {
            if (s < 1.0 - tol)
            {
                dsdT_[cellI] = slope;
                T0_[cellI] = Tsol_.value() + tol;
            }
            else
            {
                dsdT_[cellI] = 0.0;
                T0_[cellI] = T;
            }
        }
        else
        {
            dsdT_[cellI] = -1.0/(Tliq_-Tsol_).value();
            T0_[cellI] = Tsol_.value() + (s - 1.0)/dsdT_[cellI];
        }
    }

    const volScalarField::Boundary& TBoundary = T_.boundaryField();
    forAll(TBoundary, patchI)
    {
        const fvPatchField<scalar>& TPatch = TBoundary[patchI];
        const fvPatchField<scalar>& sfPatch = sf_.boundaryField()[patchI];
        fvPatchField<scalar>& T0Patch = T0_.boundaryFieldRef()[patchI];
        fvPatchField<scalar>& dsdTPatch = dsdT_.boundaryFieldRef()[patchI];
        forAll(TPatch, faceI)
        {
            const scalar T = TPatch[faceI];
            const scalar s = sfPatch[faceI];

            if (T > Tliq_.value())
            {
                if (s > tol)
                {
                    dsdTPatch[faceI] = slope;
                    T0Patch[faceI] = Tliq_.value() - tol;
                }
                else
                {
                    dsdTPatch[faceI] = 0.0;
                    T0Patch[faceI] = T;
                }
            }
            else if (T < Tsol_.value())
            {
                if (s < 1.0 - tol)
                {
                    dsdTPatch[faceI] = slope;
                    T0Patch[faceI] = Tsol_.value() + tol;
                }
                else
                {
                    dsdTPatch[faceI] = 0.0;
                    T0Patch[faceI] = T;
                }
            }
            else
            {
                dsdTPatch[faceI] = -1.0/(Tliq_-Tsol_).value();
                T0Patch[faceI] = Tsol_.value() + (s - 1.0)/dsdTPatch[faceI];
            }
        }
    }

    const dimensionedScalar rDeltaT = 1.0 / mesh_.time().deltaT();

    eqn -=
    (
        alpha_*thermo_.rho()*Lm_*fvc::ddt(sf_)
      + rDeltaT*alpha_*thermo_.rho()*Lm_*(fvm::Sp(dsdT_, T_) - dsdT_*T0_)
    );
}


// ************************************************************************* //
