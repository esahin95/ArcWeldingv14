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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::DTRMParticle::trackingData::trackingData
(
    const lagrangian::Cloud<DTRMParticle>& cloud,
    const interpolations::cellPoint<scalar>& alphaInterp,
    const interpolations::cellPoint<scalar>& absorpInterp,
    const interpolations::cellPoint<vector>& gradAlphaInterp,
    volScalarField& Q,
    DynamicField<vector>& allPositions,
    DynamicField<label>& allTracks,
    DynamicField<scalar>& allPowers,
    const meshSearch& searchEngine,
    const reflectionModel& reflection
)
:
    particle::trackingData(cloud),
    alphaInterp_(alphaInterp),
    absorpInterp_(absorpInterp),
    gradAlphaInterp_(gradAlphaInterp),
    Q_(Q),
    allPositions_(allPositions),
    allTracks_(allTracks),
    allPowers_(allPowers),
    searchEngine_(searchEngine),
    reflection_(reflection)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolations::cellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::alphaInterp() const
{
    return alphaInterp_;
}


inline const Foam::interpolations::cellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::absorpInterp() const
{
    return absorpInterp_;
}


inline const Foam::interpolations::cellPoint<Foam::vector>&
Foam::DTRMParticle::trackingData::gradAlphaInterp() const
{
    return gradAlphaInterp_;
}


inline const Foam::meshSearch&
Foam::DTRMParticle::trackingData::searchEngine() const
{
    return searchEngine_;
}


inline const Foam::reflectionModel&
Foam::DTRMParticle::trackingData::reflection() const
{
    return reflection_;
}


inline Foam::scalar& Foam::DTRMParticle::trackingData::Q(const label celli)
{
    return Q_[celli];
}


inline void Foam::DTRMParticle::trackingData::append
(
    const vector& position,
    const label trackIndex,
    const scalar power
)
{
    if (mesh.time().writeTime())
    {
        allPositions_.append(position);
        allTracks_.append(trackIndex);
        allPowers_.append(power);
    }
}


template<class Type>
inline Type Foam::DTRMParticle::interpolate
(
    const interpolations::cellPoint<Type>& interp,
    const polyMesh& mesh
) const
{
    return interp.interpolate(coordinates(), currentTetIndices(mesh));
}


// ************************************************************************* //
