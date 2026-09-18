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
    Foam::evaporationModels::Knight

Description
    Evaporation with the mass flux and recoil pressure from Knight's jump
    conditions across the Knudsen layer.

    The Mach number of the vapour increases linearly from 0 at the boiling
    temperature Tv to its maximum sqrt(gv/2) at Th. The source terms and the
    remaining parameters are those of gasDynamic.

Usage
    In addition to the gasDynamic parameters:
    \table
        Property | Description                            | Required | Default
        T0       | Ambient temperature                    | yes      |
        g0       | Heat capacity ratio of the ambient gas | yes      |
        M0       | Molar mass of the ambient gas [g/mol]  | yes      |
        Th       | Temperature of the maximum Mach number | yes      |
        gv       | Heat capacity ratio of the vapour      | yes      |
    \endtable

SourceFiles
    KnightEvaporationModel.C

See also
    Foam::evaporationModels::gasDynamic

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
                           Class Knight Declaration
\*---------------------------------------------------------------------------*/

class Knight
:
    public gasDynamic
{
protected:

    // Protected Data

        //- Ambient temperature
        const dimensionedScalar T0_;

        //- Heat capacity ratio of the ambient gas
        const dimensionedScalar g0_;

        //- Specific gas constant of the ambient gas
        const dimensionedScalar R0_;

        //- Temperature at which the Mach number reaches its maximum
        const dimensionedScalar Th_;

        //- Heat capacity ratio of the vapour
        const dimensionedScalar gv_;


public:

    //- Runtime type information
    TypeName("Knight");


    // Constructors

        //- Construct for the phase given by group
        Knight(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        Knight(const Knight&) = delete;


    //- Destructor
    virtual ~Knight()
    {}


    // Member Functions

        //- Correct the mass flux and return its maximum change
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
