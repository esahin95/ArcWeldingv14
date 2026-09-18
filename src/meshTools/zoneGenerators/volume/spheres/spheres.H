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

Class
    Foam::zoneGenerators::spheres

Description
    A zoneGenerator which selects points, cells or faces with centres either
    inside or outside of a set of spheres.

SourceFiles
    spheres.C

\*---------------------------------------------------------------------------*/

#ifndef spheres_zoneGenerator_H
#define spheres_zoneGenerator_H

#include "volume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                              Class spheres Declaration
\*---------------------------------------------------------------------------*/

class spheres
:
    public volume
{
    // Private Data

        //- Centres
        List<vector> centres_;

        //- Squared radii
        List<scalar> radiiSqr_;

        //- Translate spheres
        const vector offset_;

        //- Inflate spheres
        const scalar scale_;


    // Private Member Functions

        friend class volume;

        //- Return true if the spheres contain the given point
        inline bool contains(const point& p) const;


public:

    //- Runtime type information
    TypeName("spheres");


    // Constructors

        //- Construct from dictionary
        spheres
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~spheres();


    // Member Functions

        virtual zoneSet generate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
