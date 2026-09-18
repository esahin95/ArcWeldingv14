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

#include "evaporationModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word
Foam::evaporationModel::dictName_ = "evaporation";

namespace Foam
{
    defineTypeNameAndDebug(evaporationModel, 0);
    defineRunTimeSelectionTable(evaporationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModel::evaporationModel
(
    const fvMesh& mesh,
    const word& group
)
:
    physicalProperties(mesh, group),

    mesh_(mesh),

    thermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, group)
        )
    ),

    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", group)
        )
    ),

    rho_(mesh.lookupObject<volScalarField>("rho")),

    T_(mesh.lookupObject<volScalarField>("T")),

    alphaRhoPhi_
    (
        mesh.lookupObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", group)
        )
    ),

    mDot_
    (
        IOobject
        (
            IOobject::groupName("mDot", group),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime/dimArea, 0.0)
    )
{}


// ************************************************************************* //