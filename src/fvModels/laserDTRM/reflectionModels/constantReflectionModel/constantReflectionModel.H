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
    Foam::reflectionModels::constant

Description

SourceFiles
    constantReflectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantReflectionModel_H
#define constantReflectionModel_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reflectionModels
{


/*---------------------------------------------------------------------------*\
                          Class constantReflectionModel Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public reflectionModel
{
    // Private Data

        //- constant reflectivity
        const scalar rho_;


public:

    //- Runtime type information
    TypeName("constant")


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        constant(const constant&) = delete;


    //- Destructor
    virtual ~constant()
    {}


    // Member Functions

        // Reflectivity
        virtual scalar rho(const scalar cosTheta) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const constant&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reflectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
