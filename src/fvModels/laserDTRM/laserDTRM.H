/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::laserDTRM

Description
    Laser heating by ray tracing through a VoF interface (discrete transfer
    radiation method), derived from the implementation in OpenFOAM-v2512.

    Rays start from the positions given by the power model and are traced in
    the laser direction through the transparent phase. Where a ray enters the
    absorbing phase (alpha < 0.5) it is split: a reflected ray carrying the
    fraction given by the reflection model continues from the interface, and
    the remainder is absorbed along the ray with the absorption coefficient
    absorption*(1 - alpha). A ray is no longer traced once its power has
    fallen below 1% of its initial power. The absorbed power density is
    written as the field \c Q and added to the temperature equation; the
    traced rays are written at write times.

Usage
    \verbatim
    laser
    {
        type            laserDTRM;
        phase           gas;        // Transparent phase
        absorption      1e6;        // Absorption coefficient, default 1e6
        relax           1;          // Relaxation of Q in time, default 1
        setFormat       raw;        // Format of the ray output

        powerModel
        {
            type        random;
            position    (0 0 0);    // Beam centre, Function1 of time
            radius      1e-4;       // Beam radius
            normal      (0 -1 0);   // Beam direction
            Q           100;        // Laser power
            nRays       1000;
            powerDist   1;          // Radial power profile, Function1
        }

        reflectionModel
        {
            type        fresnelLaser;
            epsilon     0.08;
        }
    }
    \endverbatim

SourceFiles
    laserDTRM.C

\*---------------------------------------------------------------------------*/

#ifndef laserDTRM_H
#define laserDTRM_H

#include "fvModel.H"
#include "volFields.H"
#include "setWriter.H"
#include "DynamicField.H"
#include "reflectionModel.H"
#include "powerModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class laserDTRM Declaration
\*---------------------------------------------------------------------------*/

class laserDTRM
:
    public fvModel
{
    // Private Data

        //- Phase fraction of the transparent phase
        const volScalarField& alpha_;

        //- Power model generating the rays
        autoPtr<powerModel> powerModelPtr_;

        //- Time index of the last ray tracing
        label curTimeIndex_;

        //- Relaxation factor of the absorbed power between time steps
        const scalar relax_;

        //- Absorption coefficient of the absorbing phase
        const scalar a_;

        //- Reflection model
        const autoPtr<reflectionModel> reflectionModelPtr_;

        //- Absorbed power density
        volScalarField Q_;


        // Ray output

            //- Ray writer
            const autoPtr<setWriter> formatterPtr_;

            //- Output directory
            const fileName outputPath_;

            //- Positions along the rays
            mutable DynamicField<vector> allPositions_;

            //- Ray indices of the positions
            mutable DynamicField<label> allTracks_;

            //- Ray powers at the positions
            mutable DynamicField<scalar> allPowers_;

            //- Index of the last written ray file
            mutable label writeIndex_;


public:

    //- Runtime type information
    TypeName("laserDTRM");


    // Constructors

        //- Construct from explicit source name and mesh
        laserDTRM
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        laserDTRM(const laserDTRM&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source
            //  term to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Trace the rays and update the absorbed power, once per time
            //  step
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add the absorbed power to the temperature equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& T,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Write the rays traced at this time
            virtual bool write(const bool write) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laserDTRM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
