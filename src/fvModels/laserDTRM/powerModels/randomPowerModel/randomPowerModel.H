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
    Foam::powerModels::constant

Description

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
                          Class randomPowerModel Declaration
\*---------------------------------------------------------------------------*/

class random
:
    public powerModel
{
    // Private Data

        //- Random number generator
        static randomGenerator rndGen_;

        //- Number of rays
        const label nRays_;

        //- Radial power distribution
        const autoPtr<Function1<scalar>> powerDist_;


public:

    //- Runtime type information
    TypeName("random")


    // Constructors

        //- Construct from components
        random
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        random(const random&) = delete;


    //- Destructor
    virtual ~random()
    {}


    // Member Functions

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
