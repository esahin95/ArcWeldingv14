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
    Foam::solidificationModels::none

Description
    No solidification: the solid fraction remains zero and no source terms
    are added.

SourceFiles
    noSolidificationModel.C

\*---------------------------------------------------------------------------*/

#ifndef noSolidificationModel_H
#define noSolidificationModel_H

#include "solidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public solidificationModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct for the phase given by group
        none(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    virtual ~none()
    {}


    // Member Functions

        //- Correct the solid fraction; does nothing
        virtual scalar correct(const bool relax=true);

        //- Add the latent heat source; does nothing
        virtual void addSup(fvMatrix<scalar>& eqn) const;

        //- Add the solid damping source; does nothing
        virtual void addSup(fvMatrix<vector>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
