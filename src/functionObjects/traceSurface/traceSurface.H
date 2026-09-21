/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::functionObjects::traceSurface

Description
    Traces rays through a VoF phase normal to a rectangle and writes the
    layer height or fill fraction of the phase over a grid of tiles.

    The rectangle is given by three corners: the first two span its
    x-direction, the second and third its y-direction. It is divided into
    nx by ny tiles, each sampled by ns by ns rays starting on the rectangle
    and travelling along its normal x ^ y, up to the given height or the
    first wall or other non-coupled boundary. For each tile the output
    quantity is either:

      - layerHeight: the largest distance along a ray to the exit of the last
        cell with alpha >= 0.5
      - fillFraction: the mean integral of alpha along the rays, divided by
        the layer height

    At every write, one line per tile is written to the output file with the
    tile centre, in coordinates of the rectangle, and the output quantity.

Usage
    \verbatim
    layerHeight
    {
        type        traceSurface;
        libs        ("libtraceSurface.so");
        alpha       alpha.metal;
        corners     ((0 0 1e-4) (1.5e-4 0 1e-4) (1.5e-4 1e-4 1e-4));
        nx          30;
        ny          20;
        ns          8;
        height      8e-5;
        output      layerHeight;
    }
    \endverbatim

    \table
        Property | Description                            | Required | Default
        alpha    | Name of the phase fraction field       | yes      |
        corners  | Three corners of the rectangle         | yes      |
        nx       | Number of tiles in the x-direction     | yes      |
        ny       | Number of tiles in the y-direction     | yes      |
        ns       | Number of rays per tile per direction  | yes      |
        height   | Maximum ray length                     | no       | great
        output   | Output: layerHeight or fillFraction    | yes      |
    \endtable

SourceFiles
    traceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef traceSurface_H
#define traceSurface_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class traceSurface Declaration
\*---------------------------------------------------------------------------*/

class traceSurface
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Data

        //- Output quantities
        enum class outputType
        {
            layerHeight,
            fillFraction
        };

        //- Output quantity names
        static const NamedEnum<outputType, 2> outputTypeNames_;

        //- Name of the phase fraction field
        const word alphaName_;

        //- Corners of the sampled rectangle: the first two span its
        //  x-direction, the second and third its y-direction
        const List<vector> corners_;

        //- Number of tiles in the x-direction
        const label nx_;

        //- Number of tiles in the y-direction
        const label ny_;

        //- Number of rays per tile in each direction
        const label ns_;

        //- Maximum ray length
        const scalar maxTrackLength_;

        //- Output quantity
        const outputType output_;


    // Private Member Functions

        //- Trace the rays and write the output quantity of each tile
        void writePositions();


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("traceSurface");


    // Constructors

        //- Construct from Time and dictionary
        traceSurface
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~traceSurface();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute; does nothing
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
