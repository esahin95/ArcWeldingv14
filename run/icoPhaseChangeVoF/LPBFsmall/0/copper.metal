/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.copper;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [];

internalField
{
    type        zonal;

    defaultValue    0;

    zones
    {
        powderBed
        {
            type            spheres;

            data
            (
                #include "$FOAM_CASE/system/powderBed"
            );

            offset          (0 0 100e-6);
            scale           1.0;

            value           1;
        }
    }
}

boundaryField
{
    "(xmin|xmax|ymin|ymax|zmin)"
    {
        type            zeroGradient;
    }

    zmax
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
}


// ************************************************************************* //