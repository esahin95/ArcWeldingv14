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

#include "DTRMParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DTRMParticle, 0);
}

Foam::label Foam::DTRMParticle::nTracks = 0;

Foam::scalar Foam::DTRMParticle::qLost = 0;


// * * * * * * * * * * * * * * * Local Constants * * * * * * * * * * * * * * //

namespace
{
    //- Tolerance on the phase fraction when locating the interface
    const Foam::scalar alphaTol = 1e-2;

    //- Maximum number of iterations locating the interface
    const Foam::label maxInterfaceIter = 20;

    //- Fraction of the initial power below which a ray is no longer traced
    const Foam::scalar qMinFraction = 0.01;

    //- Step length below which a step is considered to have made no progress
    const Foam::scalar dsMin = 1e-10;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::DTRMParticle::locateInterface
(
    const vector& p0,
    const vector& dx,
    const label celli,
    const scalar alpha0,
    const scalar alpha1,
    const trackingData& td
)
{
    scalar a = 0, fa = alpha0 - 0.5;
    scalar b = 1, fb = alpha1 - 0.5;
    scalar c = 0, fc = fa;

    label i = 0;
    while (mag(fc) > alphaTol && i < maxInterfaceIter)
    {
        i++;

        c = (a*fb - b*fa)/(fb - fa);

        locate(td.searchEngine(), p0 + c*dx, celli);
        fc = interpolate(td.alphaInterp(), td.mesh) - 0.5;

        if (fc*fa > 0)
        {
            a = c;
            fa = fc;
        }
        else
        {
            b = c;
            fb = fc;
        }
    }

    if (i >= maxInterfaceIter)
    {
        DebugInfo<< "Failed search: " << alpha0 << " -> " << fc + 0.5 << endl;
    }

    return c;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const meshSearch& searchEngine,
    const vector& position,
    const label celli,
    label& nLocateBoundaryHits,
    const vector& direction,
    const scalar power,
    const bool transmissive
)
:
    particle(searchEngine, position, celli, nLocateBoundaryHits),
    q0_(power),
    trackIndex_(nTracks++),
    q_(power),
    d_(direction),
    transmissive_(transmissive)
{
    reset(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::DTRMParticle::move
(
    lagrangian::Cloud<DTRMParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    td.append(position(td.mesh), trackIndex_, q_);

    while
    (
        q_ > qMinFraction*q0_
     && td.keepParticle && td.sendToProc == -1 && stepFraction() < 1
    )
    {
        // State at the start of the step
        const scalar oldAlpha = interpolate(td.alphaInterp(), td.mesh);
        const scalar oldAbsorp = interpolate(td.absorpInterp(), td.mesh);
        const label oldCell = cell();
        const vector oldPos = position(td.mesh);

        if (transmissive_ && oldAlpha < 0.5 - alphaTol)
        {
            DebugInfo<< "Give up transmissive ray inside material" << endl;
            qLost += q_;
            td.keepParticle = false;
            continue;
        }

        // Track to the next face
        trackToAndHitFace(d_, 1.0, cloud, td);

        const scalar alpha = interpolate(td.alphaInterp(), td.mesh);
        const vector pos = position(td.mesh);
        const vector dx = pos - oldPos;
        const scalar ds = mag(dx);

        // Skip steps which made no progress
        if (ds < dsMin)
        {
            continue;
        }

        // On entering the absorbing phase split off the reflected ray
        if (oldAlpha >= 0.5 && alpha < 0.5 && transmissive_)
        {
            transmissive_ = false;

            label nLocateBoundaryHits = 0;
            autoPtr<DTRMParticle> reflected
            (
                new DTRMParticle
                (
                    td.searchEngine(),
                    oldPos,
                    oldCell,
                    nLocateBoundaryHits,
                    d_,
                    q0_,
                    true
                )
            );

            const scalar c = reflected->locateInterface
            (
                oldPos,
                dx,
                oldCell,
                oldAlpha,
                alpha,
                td
            );

            // Interface normal, pointing into the transparent phase
            const vector gradAlpha =
                reflected->interpolate(td.gradAlphaInterp(), td.mesh);
            const vector nHat = gradAlpha/mag(gradAlpha);

            // Split the power according to the reflectivity
            const scalar cosTheta = -nHat & d_;
            reflected->d_ = td.reflection().R(d_, nHat);
            reflected->q_ = td.reflection().rho(cosTheta)*q_;
            q_ -= reflected->q_;

            if (cosTheta < 0)
            {
                DebugInfo
                    << "Give up particle with negative angle" << nl
                    << ds << " " << oldAlpha << " " << alpha << " "
                    << (gradAlpha & d_) << " " << c << endl;

                qLost += q_ + reflected->q_;
                td.keepParticle = false;
                continue;
            }

            cloud.addParticle(reflected.ptr());
        }

        // Absorption along the step
        if (!transmissive_)
        {
            const scalar qAbsorbed = max(min(ds*oldAbsorp, 1.0), 0.0)*q_;
            td.Q(oldCell) += qAbsorbed;
            q_ -= qAbsorbed;
        }

        td.append(pos, trackIndex_, q_);
    }

    return td.keepParticle;
}


void Foam::DTRMParticle::hitWallPatch
(
    lagrangian::Cloud<DTRMParticle>&,
    trackingData& td
)
{
    td.keepParticle = false;
}


// ************************************************************************* //
