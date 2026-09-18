/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "spheres.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(spheres, 0);
        addToRunTimeSelectionTable(zoneGenerator, spheres, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::zoneGenerators::spheres::contains(const point& p) const
{
    forAll(centres_, i)
    {
        if (magSqr(centres_[i] - p) <= radiiSqr_[i])
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::spheres::spheres
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    volume(name, mesh, dict)
{
    const vector offset
    (
        dict.lookupOrDefault<vector>("offset", dimLength, Zero)
    );
    const scalar scale(dict.lookupOrDefault<scalar>("scale", dimless, 1.0));

    // Sphere data, (x y z r id) per sphere
    const label nData = 5;
    const List<scalar> data(dict.lookup("data"));

    if (data.size() % nData != 0)
    {
        FatalIOErrorInFunction(dict)
            << "size of data (" << data.size()
            << ") is not a multiple of " << nData
            << exit(FatalIOError);
    }

    centres_.setSize(data.size()/nData);
    radiiSqr_.setSize(centres_.size());

    forAll(centres_, i)
    {
        const label j = nData*i;

        centres_[i] = vector(data[j], data[j + 1], data[j + 2]) + offset;
        radiiSqr_[i] = sqr(data[j + 3]*scale);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::spheres::~spheres()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::spheres::generate() const
{
    return volume::generate(*this);
}


// ************************************************************************* //
