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
    Foam::reflectionModels::fresnelLorentzDrude

Description
    Reflectivity from the Fresnel equations for an absorbing medium, with the
    complex refractive index from the Drude free-electron model for the
    given free electron density, electrical conductivity and laser
    wavelength.

Usage
    \table
        Property | Description                   | Required | Default
        N        | Free electron density [1/m^3] | yes      |
        lambda   | Laser wavelength [m]          | yes      |
        sigma    | Electrical conductivity [S/m] | yes      |
    \endtable

SourceFiles
    fresnelLorentzDrudeReflectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef fresnelLorentzDrudeReflectionModel_H
#define fresnelLorentzDrudeReflectionModel_H

#include "reflectionModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reflectionModels
{

/*---------------------------------------------------------------------------*\
                     Class fresnelLorentzDrude Declaration
\*---------------------------------------------------------------------------*/

class fresnelLorentzDrude
:
    public reflectionModel
{
    // Private Data

        //- Free electron density
        const dimensionedScalar N_;

        //- Laser wavelength
        const dimensionedScalar lambda_;

        //- Electrical conductivity
        const dimensionedScalar sigma_;

        //- Complex refractive index, real part
        scalar n_;

        //- Complex refractive index, imaginary part
        scalar k_;


public:

    //- Runtime type information
    TypeName("fresnelLorentzDrude");


    // Constructors

        //- Construct from dictionary and mesh
        fresnelLorentzDrude(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        fresnelLorentzDrude(const fresnelLorentzDrude&) = delete;


    //- Destructor
    virtual ~fresnelLorentzDrude()
    {}


    // Member Functions

        //- Return the reflectivity for the cosine of the angle of incidence
        virtual scalar rho(const scalar cosTheta) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fresnelLorentzDrude&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reflectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
