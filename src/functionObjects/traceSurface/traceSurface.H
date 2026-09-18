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
    // Private enumeration

        //- Output enumeration
        enum class outputType
        {
            layerHeight,
            fillFraction
        };

        //- Output enumeration names
        static const NamedEnum<outputType, 2> outputTypeNames_;

    // Private Data

        //- Name of the alpha field
        const word alphaName_;

        //- Area to sample
        const List<vector> corners_;

        //- Discretization
        const label nx_;
        const label ny_;
        const label ns_;

        //- Track length
        const scalar maxTrackLength_;

        //- Output type
        const outputType output_;


    // Private Member Functions

        //- Output positions
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

        //- Execute
        virtual bool execute();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
