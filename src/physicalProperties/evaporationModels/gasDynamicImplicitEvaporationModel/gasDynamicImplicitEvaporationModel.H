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
    Foam::gasDynamicImplicitEvaporationModel

Description

SourceFiles
    gasDynamicImplicitEvaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef gasDynamicImplicitEvaporationModel_H
#define gasDynamicImplicitEvaporationModel_H

#include "gasDynamicEvaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace evaporationModels
{


/*---------------------------------------------------------------------------*\
                          Class gasDynamicImplicitEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class gasDynamicImplicit
:
    public gasDynamic
{

public:

    //- Runtime type information
    TypeName("gasDynamicImplicit");

    // Constructors

        //- Construct from components
        gasDynamicImplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        gasDynamicImplicit(const gasDynamicImplicit&) = delete;


    //- Destructor
    virtual ~gasDynamicImplicit()
    {}


    // Member functions

        //- Source term for energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gasDynamicImplicit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
