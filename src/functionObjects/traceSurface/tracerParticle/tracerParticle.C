/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "tracerParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tracerParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tracerParticle::tracerParticle
(
    const meshSearch& searchEngine,
    const vector& position,
    const label celli,
    label& nLocateBoundaryHits,
    const vector& direction
)
:
    particle(searchEngine, position, celli, nLocateBoundaryHits),
    h_(0),
    d_(direction),
    a_(0),
    s_(0)
{
    reset(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::tracerParticle::move
(
    lagrangian::Cloud<tracerParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    while (td.keepParticle && td.sendToProc == -1 && stepFraction() < 1)
    {
        const vector pos0 = position(td.mesh);
        const scalar alpha = td.alpha(cell());

        // Track to the next face
        trackToAndHitFace(d_, 1.0, cloud, td);

        const scalar ds = mag(position(td.mesh) - pos0);

        s_ += ds;
        a_ += alpha*ds;

        if (alpha >= 0.5)
        {
            h_ = s_;
        }
    }

    return td.keepParticle;
}


void Foam::tracerParticle::hitBasicPatch
(
    lagrangian::Cloud<tracerParticle>&,
    trackingData&
)
{
    stepFraction() = 1;
}


void Foam::tracerParticle::hitWallPatch
(
    lagrangian::Cloud<tracerParticle>&,
    trackingData&
)
{
    stepFraction() = 1;
}


// ************************************************************************* //
