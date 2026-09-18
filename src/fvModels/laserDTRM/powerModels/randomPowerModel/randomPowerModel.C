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

#include "randomPowerModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace powerModels
    {
        defineTypeNameAndDebug(random, 0);
        addToRunTimeSelectionTable
        (
            powerModel,
            random,
            dictionary
        );

        randomGenerator random::rndGen_(261782, true);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::random::random
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    powerModel(dict, mesh),

    nRays_(dict.lookup<scalar>("nRays")),

    powerDist_
    (
        Function1<scalar>::New
        (
            "powerDist",
            dimless,
            dimless,
            dict
        )
    )
{
    // Resize lists
    positions_.resize(nRays_);
    powers_.resize(nRays_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::powerModels::random::initialise()
{
    // Current laser position
    const vector centre = pos_->value(mesh_.time().value());
    //const vector centre = pos_->value(0.0);

    // Right handed CS
    const vector t1 = normalised(perpendicular(normal_));
    const vector t2 = normalised(normal_ ^ t1);

    forAll(positions_, i)
    {
        // Rejection sampling from arbitrary distribution
        vector delta;
        scalar deltaMag;
        do
        {
            delta =
                  t1*rndGen_.scalarAB(-rad_, rad_)
                + t2*rndGen_.scalarAB(-rad_, rad_);

            deltaMag = mag(delta);
        }
        while
        (
            deltaMag > rad_ ||
            rndGen_.scalar01() > powerDist_->value(deltaMag)
        );
        positions_[i] = centre + delta;

        // Constant power per ray
        powers_[i] = Q_ / nRays_;
    }
}


// ************************************************************************* //
