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

#include "constantReflectionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace reflectionModels
    {
        defineTypeNameAndDebug(constant, 0);
        addToRunTimeSelectionTable
        (
            reflectionModel,
            constant,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectionModels::constant::constant
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    reflectionModel(dict, mesh),
    rho_(dict.lookup<scalar>("rho"))
{
    DebugInfo<< "Reflectivity set to rho = " << rho_ <<endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::reflectionModels::constant::rho
(
    const scalar cosTheta
) const
{
    return rho_;
}


// ************************************************************************* //
