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
    Foam::KnightEvaporationModel

Description

SourceFiles
    KnightEvaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef KnightEvaporationModel_H
#define KnightEvaporationModel_H

#include "gasDynamicEvaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace evaporationModels
{


/*---------------------------------------------------------------------------*\
                          Class KnightEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class Knight
:
    public gasDynamic
{


protected:

    // Evaporation properties

        //- Atmospheric properties
        const dimensionedScalar T0_;
        const dimensionedScalar g0_;
        const dimensionedScalar R0_;

        //- Vapor properties
        const dimensionedScalar Th_;
        const dimensionedScalar gv_;


public:

    //- Runtime type information
    TypeName("Knight");

    // Constructors

        //- Construct from components
        Knight(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        Knight(const Knight&) = delete;


    //- Destructor
    virtual ~Knight()
    {}


    // Member functions

        //- Correct solid fraction
        virtual scalar correct(const bool relax=true);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Knight&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
