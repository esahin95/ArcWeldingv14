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

Class
    Foam::powerModels::random

Description
    Rays starting at random positions within the beam radius, each carrying
    the power Q/nRays.

    The positions are sampled by rejection with the acceptance probability
    powerDist(r), a Function1 of the distance r from the beam centre, so the
    power per unit area is proportional to powerDist. The same random
    sequence is used on all processors.

Usage
    \table
        Property  | Description                       | Required | Default
        position  | Beam centre, Function1 of time    | yes      |
        radius    | Beam radius                       | yes      |
        normal    | Beam direction                    | yes      |
        Q         | Laser power                       | yes      |
        nRays     | Number of rays                    | yes      |
        powerDist | Radial power profile, Function1   | yes      |
    \endtable

SourceFiles
    randomPowerModel.C

\*---------------------------------------------------------------------------*/

#ifndef randomPowerModel_H
#define randomPowerModel_H

#include "powerModel.H"
#include "randomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class random Declaration
\*---------------------------------------------------------------------------*/

class random
:
    public powerModel
{
    // Private Data

        //- Random number generator, synchronised across processors
        static randomGenerator rndGen_;

        //- Number of rays
        const label nRays_;

        //- Radial power profile
        const autoPtr<Function1<scalar>> powerDist_;


public:

    //- Runtime type information
    TypeName("random");


    // Constructors

        //- Construct from dictionary and mesh
        random(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        random(const random&) = delete;


    //- Destructor
    virtual ~random()
    {}


    // Member Functions

        //- Sample new start positions for the current beam position
        virtual void initialise();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const random&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
