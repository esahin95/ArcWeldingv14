/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      binary;
    class       volScalarField;
    location    "0";
    object      T.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 300;

boundaryField
{
    ymax
    {
        type            calculated;
        value           uniform 300;
    }
    xmin
    {
        type            calculated;
        value           uniform 300;
    }
    xmax
    {
        type            calculated;
        value           uniform 300;
    }
    ymin
    {
        type            calculated;
        value           uniform 300;
    }
    zmin
    {
        type            calculated;
        value           uniform 300;
    }
    zmax
    {
        type            calculated;
        value           uniform 300;
    }
}


// ************************************************************************* //
