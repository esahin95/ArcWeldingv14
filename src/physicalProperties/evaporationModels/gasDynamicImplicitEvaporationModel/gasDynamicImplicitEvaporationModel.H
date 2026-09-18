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
    Foam::evaporationModels::gasDynamicImplicit

Description
    Evaporation in the gas-dynamic limit as gasDynamic, with the evaporative
    cooling linearised in the temperature.

    Using d(mDot)/dT = mDot/T*(L/(Rv*T) - 1/2), the cooling L*mDot*|grad(alpha)|
    is split into a part implicit in T and an explicit remainder. The recoil
    pressure and the parameters are those of gasDynamic.

SourceFiles
    gasDynamicImplicitEvaporationModel.C

See also
    Foam::evaporationModels::gasDynamic

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
                     Class gasDynamicImplicit Declaration
\*---------------------------------------------------------------------------*/

class gasDynamicImplicit
:
    public gasDynamic
{
public:

    //- Runtime type information
    TypeName("gasDynamicImplicit");


    // Constructors

        //- Construct for the phase given by group
        gasDynamicImplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        gasDynamicImplicit(const gasDynamicImplicit&) = delete;


    //- Destructor
    virtual ~gasDynamicImplicit()
    {}


    // Member Functions

        //- Recoil pressure as gasDynamic
        using gasDynamic::addSup;

        //- Add the linearised evaporative cooling to the energy equation
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
