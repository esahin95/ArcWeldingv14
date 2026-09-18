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
    Must be updated. --- IGNORE ---
    This function object reports the height of the interface above a set of
    locations. For each location, it writes the vertical distance of the
    interface above both the location and the lowest boundary. It also writes
    the point on the interface from which these heights are computed. It uses
    an integral approach, so if there are multiple interfaces above or below a
    location then this method will generate average values.

    Example of function object specification:
    \verbatim
    traceSurface1
    {
        type           traceSurface;
        libs           ("libfieldFunctionObjects.so");
        alpha          alpha.water;
        locations      ((0 0 0) (10 0 0) (20 0 0));
    }
    \endverbatim

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
