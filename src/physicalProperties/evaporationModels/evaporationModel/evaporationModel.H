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
    Foam::evaporationModel

Description
    Base class for evaporation models of a VoF phase.

    The model is selected from the optional \c evaporation sub-dictionary of
    the phase physicalProperties; without it no evaporation is modelled
    (\c none). A model maintains the evaporation mass flux per unit
    interface area, \c mDot.<phase>, and provides the corresponding source
    terms for the energy and momentum equations.

Usage
    In constant/physicalProperties.<phase>:
    \verbatim
    evaporation
    {
        type    gasDynamic;
        L       6.3e6;
        p0      1e5;
        Mv      56;
        Tv      3100;
        relax   0.9;
    }
    \endverbatim

SourceFiles
    evaporationModel.C
    evaporationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef evaporationModel_H
#define evaporationModel_H

#include "physicalProperties.H"
#include "volFields.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class evaporationModel
:
    public physicalProperties
{
protected:

    // Protected Data

        //- Name of the model sub-dictionary
        static const word dictName_;

        //- Phase fraction
        const volScalarField& alpha_;

        //- Temperature
        const volScalarField& T_;

        //- Evaporation mass flux per unit interface area
        volScalarField::Internal mDot_;


public:

    //- Runtime type information
    TypeName("evaporationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            evaporationModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& group
            ),
            (mesh, group)
        );


    // Constructors

        //- Construct for the phase given by group
        evaporationModel(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        evaporationModel(const evaporationModel&) = delete;


    // Selectors

        //- Select the evaporation model of the phase given by group
        static autoPtr<evaporationModel> New
        (
            const fvMesh& mesh,
            const word& group = word::null
        );


    //- Destructor
    virtual ~evaporationModel()
    {}


    // Member Functions

        //- Correct the mass flux and return its maximum change
        virtual scalar correct(const bool relax=true) = 0;

        //- Add the evaporative cooling to the energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const = 0;

        //- Add the recoil pressure to the momentum equation
        virtual void addSup(fvMatrix<vector>& eqn) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const evaporationModel&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
