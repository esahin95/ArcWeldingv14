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
    Foam::tracerParticle

Description
    Laser particle class with coupling to the phase field.

SourceFiles
    tracerParticleI.H
    tracerParticle.C
    tracerParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef tracerParticle_H
#define tracerParticle_H

#include "particle.H"
#include "Cloud.H"
#include "IOstream.H"
//#include "autoPtr.H"
//#include "interpolationCellPoint.H"
#include "volFields.H"
#include "contiguous.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class tracerParticle;

Ostream& operator<<(Ostream&, const tracerParticle&);


/*---------------------------------------------------------------------------*\
                        Class tracerParticle Declaration
\*---------------------------------------------------------------------------*/

class tracerParticle
:
    public particle
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Layer height
        scalar h_;

        //- Direction
        vector d_;

        //- Filled fraction
        scalar a_;


public:

    friend class lagrangian::Cloud<tracerParticle>;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const volScalarField& alpha_;

    public:

        // Constructors

            inline trackingData
            (
                const lagrangian::Cloud<tracerParticle>& cloud,
                const volScalarField& alpha
            );


        // Member Functions

            inline scalar alpha(label celli) const;
    };


    //- Runtime type information
    TypeName("tracerParticle");


    // Constructors

        //- Construct from components
        tracerParticle
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label cellI,
            label& nLocateBoundaryHits,
            const vector& direction
        );

        //- Construct from Istream
        tracerParticle(Istream& is, bool readFields = true);

        //- Construct from Istream and return
        static autoPtr<tracerParticle> New(Istream& is)
        {
            return autoPtr<tracerParticle>(new tracerParticle(is));
        }


    // Member Functions


        // Tracking

            //- Move
            bool move(lagrangian::Cloud<tracerParticle>&, trackingData&);

            //- Height
            scalar h() const
            {
                return h_;
            }

            //- Fill fraction
            scalar a() const
            {
                return a_;
            }

        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitBasicPatch(lagrangian::Cloud<tracerParticle>& cloud, trackingData& td);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const tracerParticle&);
};


template<>
inline bool contiguous<tracerParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tracerParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
