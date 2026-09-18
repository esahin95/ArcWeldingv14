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
    Foam::evaporationModels::gasDynamic

Description
    Evaporation in the gas-dynamic (sonic) limit.

    The saturation pressure follows the Clausius-Clapeyron relation

        pSat = p0*exp(L/(Rv*Tv)*(1 - Tv/T))

    and the evaporation mass flux and recoil pressure are

        mDot = 0.816*pSat/sqrt(2*pi*Rv*T),  pRec = 0.54*pSat

    The evaporative cooling L*mDot*|grad(alpha)| is added explicitly to the
    energy equation and the recoil pressure pRec*grad(alpha) to the momentum
    equation.

Usage
    \table
        Property | Description                          | Required | Default
        L        | Latent heat of evaporation           | yes      |
        p0       | Ambient pressure                     | yes      |
        Mv       | Molar mass of the vapour [g/mol]     | yes      |
        Tv       | Boiling temperature at p0            | yes      |
        relax    | Mass flux relaxation factor          | yes      |
    \endtable

SourceFiles
    gasDynamicEvaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef gasDynamicEvaporationModel_H
#define gasDynamicEvaporationModel_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{

/*---------------------------------------------------------------------------*\
                         Class gasDynamic Declaration
\*---------------------------------------------------------------------------*/

class gasDynamic
:
    public evaporationModel
{
protected:

    // Protected Data

        //- Model coefficients dictionary
        const dictionary& dict_;

        //- Latent heat of evaporation
        const dimensionedScalar Lv_;

        //- Ambient pressure
        const dimensionedScalar p0_;

        //- Specific gas constant of the vapour
        const dimensionedScalar Rv_;

        //- Boiling temperature at the ambient pressure
        const dimensionedScalar Tv_;

        //- Mass flux relaxation factor
        const scalar relax_;

        //- Recoil pressure
        volScalarField::Internal pRec_;


public:

    //- Runtime type information
    TypeName("gasDynamic");


    // Constructors

        //- Construct for the phase given by group
        gasDynamic(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        gasDynamic(const gasDynamic&) = delete;


    //- Destructor
    virtual ~gasDynamic()
    {}


    // Member Functions

        //- Correct the mass flux and return its maximum change
        virtual scalar correct(const bool relax=true);

        //- Add the evaporative cooling to the energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const;

        //- Add the recoil pressure to the momentum equation
        virtual void addSup(fvMatrix<vector>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gasDynamic&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
