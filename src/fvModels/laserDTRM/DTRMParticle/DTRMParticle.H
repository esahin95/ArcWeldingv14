/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::DTRMParticle

Description
    Ray of the laserDTRM model, tracked through the mesh as a particle.

    A transmissive ray crosses the transparent phase without losing power.
    On entering the absorbing phase it becomes absorbing and splits off a
    reflected, transmissive ray starting on the interface. An absorbing ray
    deposits its power in the cells it crosses.

SourceFiles
    DTRMParticleI.H
    DTRMParticle.C
    DTRMParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef DTRMParticle_H
#define DTRMParticle_H

#include "particle.H"
#include "Cloud.H"
#include "cellPoint_interpolation.H"
#include "contiguous.H"
#include "meshSearch.H"
#include "DynamicField.H"
#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class DTRMParticle;
Ostream& operator<<(Ostream&, const DTRMParticle&);

/*---------------------------------------------------------------------------*\
                        Class DTRMParticle Declaration
\*---------------------------------------------------------------------------*/

class DTRMParticle
:
    public particle
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Initial power of the ray
        scalar q0_;

        //- Index of the ray, for output
        label trackIndex_;

        //- Current power
        scalar q_;

        //- Direction
        vector d_;

        //- Is the ray transmissive, i.e. not yet absorbing?
        bool transmissive_;


public:

    //- Number of rays generated in the current ray tracing
    static label nTracks;

    //- Power lost during the ray tracing
    static scalar qLost;

    friend class lagrangian::Cloud<DTRMParticle>;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for the continuous phase fields

            //- Phase fraction of the transparent phase
            const interpolations::cellPoint<scalar>& alphaInterp_;

            //- Absorption coefficient
            const interpolations::cellPoint<scalar>& absorpInterp_;

            //- Phase fraction gradient
            const interpolations::cellPoint<vector>& gradAlphaInterp_;

        //- Absorbed power per cell
        volScalarField& Q_;

        // Ray output

            //- Positions along the rays
            DynamicField<vector>& allPositions_;

            //- Ray indices of the positions
            DynamicField<label>& allTracks_;

            //- Ray powers at the positions
            DynamicField<scalar>& allPowers_;

        //- Mesh search engine
        const meshSearch& searchEngine_;

        //- Reflection model
        const reflectionModel& reflection_;


    public:

        // Constructors

            inline trackingData
            (
                const lagrangian::Cloud<DTRMParticle>& cloud,
                const interpolations::cellPoint<scalar>& alphaInterp,
                const interpolations::cellPoint<scalar>& absorpInterp,
                const interpolations::cellPoint<vector>& gradAlphaInterp,
                volScalarField& Q,
                DynamicField<vector>& allPositions,
                DynamicField<label>& allTracks,
                DynamicField<scalar>& allPowers,
                const meshSearch& searchEngine,
                const reflectionModel& reflection
            );


        // Member Functions

            inline const interpolations::cellPoint<scalar>&
                alphaInterp() const;

            inline const interpolations::cellPoint<scalar>&
                absorpInterp() const;

            inline const interpolations::cellPoint<vector>&
                gradAlphaInterp() const;

            inline const meshSearch& searchEngine() const;

            inline const reflectionModel& reflection() const;

            //- Return the absorbed power of the given cell
            inline scalar& Q(const label celli);

            //- Record a ray position at write times
            inline void append
            (
                const vector& position,
                const label trackIndex,
                const scalar power
            );
    };


private:

    // Private Member Functions

        //- Interpolate a field to the current position
        template<class Type>
        inline Type interpolate
        (
            const interpolations::cellPoint<Type>& interp,
            const polyMesh& mesh
        ) const;

        //- Move to the interface, alpha = 0.5, on the segment p0 to p0 + dx,
        //  given the phase fractions alpha0 and alpha1 at its ends, by
        //  regula falsi. Return the fraction of dx at which it was found.
        scalar locateInterface
        (
            const vector& p0,
            const vector& dx,
            const label celli,
            const scalar alpha0,
            const scalar alpha1,
            const trackingData& td
        );


public:

    //- Runtime type information
    TypeName("DTRMParticle");


    // Constructors

        //- Construct from position and cell, direction and power
        DTRMParticle
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits,
            const vector& direction,
            const scalar power,
            const bool transmissive = true
        );

        //- Construct from Istream
        DTRMParticle(Istream& is, bool readFields = true);

        //- Construct from Istream and return
        static autoPtr<DTRMParticle> New(Istream& is)
        {
            return autoPtr<DTRMParticle>(new DTRMParticle(is));
        }


    // Member Functions

        // Tracking

            //- Track the ray, splitting off the reflected ray at the
            //  interface and depositing the absorbed power
            bool move(lagrangian::Cloud<DTRMParticle>&, trackingData&);


        // Patch interactions

            //- Stop the ray on hitting a wall
            void hitWallPatch
            (
                lagrangian::Cloud<DTRMParticle>& cloud,
                trackingData& td
            );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const DTRMParticle&);
};


template<>
inline bool contiguous<DTRMParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DTRMParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
