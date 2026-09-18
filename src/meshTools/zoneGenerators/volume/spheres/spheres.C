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
    bool result = false;
    forAll(centres_, sphereI)
    {
        result =
            result || (magSqr(centres_[sphereI] - p) <= radiiSqr_[sphereI]);

        if (result)
        {
            break;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::spheres::spheres
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    volume(name, mesh, dict),
    centres_(),
    radiiSqr_(),
    offset_(dict.lookupOrDefault<vector>("offset", dimLength, Zero)),
    scale_(dict.lookupOrDefault<scalar>("scale", dimless, 1.0))
{
    // Raw list of data (x y z r clumpID)
    List<scalar> data(dict.lookup("data"));
    if (data.size() % 5 != 0)
    {
        FatalIOErrorInFunction(dict)
            << "size of data (" << data.size()
            << ") is not a non-zero multiple of 5"
            << exit(FatalIOError);
    }

    // Read in spheres
    for(label i = 0; i < data.size(); i += 5)
    {
        centres_.append
            (
                vector(data[i], data[i+1], data[i+2]) + offset_
            );
        radiiSqr_.append(sqr(data[i+3] * scale_));
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
