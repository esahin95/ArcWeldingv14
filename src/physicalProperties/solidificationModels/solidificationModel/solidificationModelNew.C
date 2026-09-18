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

\*---------------------------------------------------------------------------*/

#include "solidificationModel.H"
#include "noSolidificationModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidificationModel> Foam::solidificationModel::New
(
    const fvMesh& mesh,
    const word& group
)
{
    const IOdictionary dict
    (
        solidificationModel::findModelDict(mesh, group)
    );

    if (dict.isDict(dictName_))
    {
        const dictionary& modelDict = dict.subDict(dictName_);

        const word modelType(modelDict.lookup("type"));

        Info<< "Selecting solidification model " << modelType << endl;

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown solidification model " << modelType << nl << nl
                << "Valid solidification models are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<solidificationModel>(cstrIter()(mesh, group));
    }
    else
    {
        Info<<"There is no " << dictName_ << " dictionary"<<endl;
        Info<<"Selecting default solidification model none"<<endl;
        return autoPtr<solidificationModel>
        (
            new solidificationModels::none(mesh, group)
        );
    }
}


// ************************************************************************* //