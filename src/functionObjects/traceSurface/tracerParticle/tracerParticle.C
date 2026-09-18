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
    const label cellI,
    label& nLocateBoundaryHits,
    const vector& direction
)
:
    particle(searchEngine, position, cellI, nLocateBoundaryHits),
    h_(0.0),
    d_(direction),
    a_(0.0)
{
    reset(0.0);
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

    scalar h = h_;

    while
    (
        td.keepParticle && td.sendToProc == -1 && stepFraction() < 1
    )
    {
        // Start position
        const vector pos0 = this->position(td.mesh);
        const label celli = cell();
        const scalar alpha = td.alpha(celli);

        // Track to new face and cell
        trackToAndHitFace(d_, 1.0, cloud, td);
        const vector pos = this->position(td.mesh);
        const scalar ds = mag(pos - pos0);

        h += ds;
        a_ += alpha*ds;

        // Update max height
        if (alpha >= 0.5)
        {
            h_ = h;
        }
    }

    return td.keepParticle;
}

void Foam::tracerParticle::hitBasicPatch
(
    lagrangian::Cloud<tracerParticle>& cloud,
    trackingData& td
)
{
    //td.keepParticle = false;
    stepFraction() = 1.0;
}


// ************************************************************************* //