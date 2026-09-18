/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight   55.8; // Predominantly iron
    }

    equationOfState
    {
        rho         8000.0;
    }

    thermodynamics
    {
        Cp          520;
        hf          0;
    }

    transport
    {
        mu          0.004; // nu = 5e-7
        kappa       10.0;
    }
}

evaporation
{
    type    gasDynamicImplicit;
    L       7.45e6;
    p0      1e5;
    Tv      3068.0;
    Mv      55.8;
    relax   1.0;//0.6;

    Th      4210.0;
    T0      300;
    M0      40.0;
    g0      1.6667;
    gv      1.6667;
}


// ************************************************************************* //
