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

#include "fresnelLaserReflectionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace reflectionModels
    {
        defineTypeNameAndDebug(fresnelLaser, 0);
        addToRunTimeSelectionTable
        (
            reflectionModel,
            fresnelLaser,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectionModels::fresnelLaser::fresnelLaser
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    reflectionModel(dict, mesh),
    epsilon_(dict.lookup<scalar>("epsilon"))
{
    DebugInfo<< "Model constant set to epsilon = " << epsilon_ <<endl;
    DebugInfo<< "Normal incidence reflectivity R = " << rho(1.0) << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::reflectionModels::fresnelLaser::rho
(
    const scalar cosTheta
) const
{
    return 0.5 *
        (
            (1 + Foam::sqr(1 - epsilon_ * cosTheta))
          / (1 + Foam::sqr(1 + epsilon_ * cosTheta))
          + (Foam::sqr(epsilon_ - cosTheta) + Foam::sqr(cosTheta))
          / (Foam::sqr(epsilon_ + cosTheta) + Foam::sqr(cosTheta))
        );
}


// ************************************************************************* //
