/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "traceSurface.H"
#include "tracerParticle.H"
#include "meshSearch.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(traceSurface, 0);
    addToRunTimeSelectionTable(functionObject, traceSurface, dictionary);
}
}

const Foam::NamedEnum<Foam::functionObjects::traceSurface::outputType, 2>
Foam::functionObjects::traceSurface::outputTypeNames_
{
    "layerHeight",
    "fillFraction"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::traceSurface::writePositions()
{
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    // Axes of the sampled rectangle and the ray direction normal to it
    const vector p0 = corners_[0];
    const vector p1 = corners_[1] - corners_[0];
    const vector p2 = corners_[2] - corners_[1];
    const vector d1 = normalised(p1);
    const vector d2 = normalised(p2);
    const vector d3 = normalised(d1 ^ d2);

    // Tile size and ray spacing within a tile
    const scalar dx = mag(p1)/scalar(nx_);
    const scalar dy = mag(p2)/scalar(ny_);
    const scalar dlx = dx/scalar(ns_);
    const scalar dly = dy/scalar(ns_);

    // Ray start points in the first tile
    List<vector> locations(ns_*ns_, Zero);
    {
        label idx = 0;
        scalar x = 0.5*dlx;
        for (label i=0; i<ns_; i++)
        {
            scalar y = 0.5*dly;
            for (label j=0; j<ns_; j++)
            {
                locations[idx++] = p0 + x*d1 + y*d2;
                y += dly;
            }
            x += dlx;
        }
    }

    const meshSearch& searchEngine = meshSearch::New(mesh());

    lagrangian::Cloud<tracerParticle> cloud
    (
        mesh(),
        "TracerCloud",
        IDLList<tracerParticle>()
    );

    tracerParticle::trackingData td(cloud, alpha);

    label nLocateBoundaryHits = 0;

    scalar x = 0;
    for (label i=0; i<nx_; i++)
    {
        scalar y = 0;
        for (label j=0; j<ny_; j++)
        {
            // Start the rays of this tile on the processors containing them
            cloud.clear();

            const vector pc = x*d1 + y*d2;
            forAll(locations, idx)
            {
                const vector p = locations[idx] + pc;
                const label celli = searchEngine.findCell(p);

                if (celli != -1)
                {
                    cloud.addParticle
                    (
                        new tracerParticle
                        (
                            searchEngine,
                            p,
                            celli,
                            nLocateBoundaryHits,
                            d3*maxTrackLength_
                        )
                    );
                }
            }

            cloud.move(cloud, td);

            // Layer height of the tile, and the sum of the filled lengths
            scalar h = 0;
            scalar sumA = 0;
            forAllConstIter(lagrangian::Cloud<tracerParticle>, cloud, iter)
            {
                h = max(h, iter().h());
                sumA += iter().a();
            }
            reduce(h, maxOp<scalar>());

            scalar out = h;

            // Fill fraction: mean filled length relative to the layer height
            if (output_ == outputType::fillFraction)
            {
                const label nRays = returnReduce(cloud.size(), sumOp<label>());
                reduce(sumA, sumOp<scalar>());

                out = h > 0 ? sumA/scalar(nRays)/h : 0;
            }

            if (Pstream::master())
            {
                const Omanip<int> w = valueWidth(1);

                file() << w << x + 0.5*dx << w << y + 0.5*dy << w << out;
                file().endl();
            }

            y += dy;
        }

        x += dx;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::traceSurface::writeFileHeader(const label i)
{
    writeHeaderValue(file(), "traced surface for ", alphaName_);

    const Omanip<int> w = valueWidth(1);
    file() << w << "# x" << w << "y" << w << "out";
    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::traceSurface::traceSurface
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(mesh_, name),
    alphaName_(dict.lookup("alpha")),
    corners_(dict.lookup("corners")),
    nx_(dict.lookup<label>("nx")),
    ny_(dict.lookup<label>("ny")),
    ns_(dict.lookup<label>("ns")),
    maxTrackLength_(dict.lookupOrDefault<scalar>("height", great)),
    output_(outputTypeNames_[dict.lookup<word>("output")])
{
    read(dict);

    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::traceSurface::~traceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::traceSurface::fields() const
{
    return wordList(alphaName_);
}


bool Foam::functionObjects::traceSurface::execute()
{
    return true;
}


bool Foam::functionObjects::traceSurface::write()
{
    logFiles::write();

    writePositions();

    return true;
}


// ************************************************************************* //
