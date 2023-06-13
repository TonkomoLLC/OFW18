    const fvMesh& mesh = fluidRegions[i];

    combustionModel& reaction = reactionFluid[i];

    fluidReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::momentumTransportModel& turbulence = turbulenceFluid[i];
    fluidReactionThermophysicalTransportModel& thermophysicalTransport =
        thermophysicalTransportFluid[i];

    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const dimensionedScalar& pRef = pRefFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    Foam::fvModels& fvModels = fvModelsFluid[i];
    Foam::fvConstraints& fvConstraints = fvConstraintsFluid[i];

    #include "checkRadiationModel.H"

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimple = pimples.pimple(i);

    pressureReference& pressureReference = pressureReferenceFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];

    // This solver does not support moving mesh but it uses the pressure
    // equation of one which does, so we need a dummy face-momentum field
    autoPtr<surfaceVectorField> rhoUf(nullptr);

// Overall U for 1D calculation

// Setup porous properties
    IOdictionary porousPropertiesDict
    (
        IOobject
        (
            "porousProperties",
            thermo.p().mesh().time().constant(),
            thermo.p().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // dimensions
    const dimensionedScalar UnitDensity("UnitDensity", dimMass/dimLength/dimLength/dimLength, 1.0);
    const dimensionedScalar UnitTemperature("UnitTemperature", dimTemperature, 1.0);
    const dimensionedScalar UnitHXCoeff("UnitHXCoeff", dimEnergy/dimTime/dimArea/dimTemperature, 1.0);
    const dimensionedScalar UnitLength("UnitLength", dimLength, 1.0);
    const dimensionedScalar UnitMeterPerSecond("UnitMeterPerSecond", dimLength/dimTime, 1.0);

    volScalarField& overallU = overallUFluid[i];
    volScalarField& kappaEffSolid = kappaEffSolidFluid[i];
    volScalarField& kappaEff = kappaEffFluid[i];
    volScalarField& alphaEff = alphaEffFluid[i];
    volScalarField& Bi = BiTubeFluid[i];    
    volScalarField& OneOverU = OneOverUFluid[i];    
    overallU *= 0;
     
    // solid density, kg/m3
    scalar rhoSolid = porousPropertiesDict.lookupOrDefault<scalar>("rhoSolid",2000);

    // solid heat capacity J/kg/K
    scalar cPSolid = porousPropertiesDict.lookupOrDefault<scalar>("cPSolid",900);

    // solid thermal conductivity W/m/K
    scalar kappaSolid = porousPropertiesDict.lookupOrDefault<scalar>("kappaSolid",5);

    // porosity
    scalar porosity = porousPropertiesDict.lookupOrDefault<scalar>("porosity",0.4);

    // porous zone name
    word porousZoneName = porousPropertiesDict.lookupOrDefault<word>("porousZoneName","porosity");    

    // Switch - use overall heat transfer coefficient
    bool useOverallU = porousPropertiesDict.lookupOrDefault<Switch>("useOverallU", false);

    // Switch - set constant overall heat transfer coefficient
    bool useOverallUconst = porousPropertiesDict.lookupOrDefault<Switch>("useOverallUconst", false);

    // Switch - if useOverallUconst is true, then overallUconst is used as a fixed overall heat transfer coefficient value, W/m2/K
    scalar overallUconst = porousPropertiesDict.lookupOrDefault<scalar>("overallUconst",100);

    // Switch - if useUSuperficial is true, then use a constant value for superficial velocity (rather than calculating mag(U)
    bool useUSuperficial = porousPropertiesDict.lookupOrDefault<Switch>("useUSuperficial", false);

    // Switch - if useUSuperficial is true, use USuperficial, m/s, for overall heat transfer coefficient calculations
    scalar USuperficial = porousPropertiesDict.lookupOrDefault<scalar>("USuperficial",1);
    
    // Switch - if use porousThermalConductivity
    bool porousThermalConductivity = porousPropertiesDict.lookupOrDefault<Switch>("porousThermalConductivity", false);
    if (useOverallUconst)  // if overallU is calculated, then porousThermalConductivity is automatically on
        porousThermalConductivity = true;
    
    // tube diameter, m
    scalar diamTube = porousPropertiesDict.lookupOrDefault<scalar>("diamTube", 0.04);

    // particle diameter, m
    scalar diamParticle = porousPropertiesDict.lookupOrDefault<scalar>("diamParticle",0.001);

    // wall temperature, K
    scalar Twall = porousPropertiesDict.lookupOrDefault<scalar>("Twall",500);

    // average density, kg/m3
    volScalarField rhoAvg = UnitDensity*rhoSolid*(1-porosity) + rho*porosity;

    // calculation of effective heat capacity
    volScalarField& CpEff = CpEffFluid[i];
    volScalarField cPFluid = reactionFluid[i].thermo().Cp();

    labelList cellsZone;

    if (useOverallU || porousThermalConductivity)
    {    
        label zoneID = mesh.cellZones().findZoneID(porousZoneName);
        const cellZone& cellsPorosity = mesh.cellZones()[zoneID];
        const meshCellZones& zoneMesh = cellsPorosity.meshZones();
        cellsZone = zoneMesh[zoneID];   //list of all porosityZone cell ID's

        forAll(rho, celli)
        {
            forAll(cellsZone, cI)
            {
                if (cellsZone[cI] == celli)
                {    
                    CpEff[celli] = ((cPSolid*(1-porosity)*rhoSolid + cPFluid[celli]*porosity*rho[celli])/rhoAvg[celli]);
                }
            }
        }

        forAll(rho.boundaryField(), patchi)
        {
            forAll(rho.boundaryField()[patchi], facei)
            {
                CpEff.boundaryFieldRef()[patchi][facei] = ((cPSolid*(1-porosity)*rhoSolid 
                    + cPFluid.boundaryField()[patchi][facei]*porosity*
                        rho.boundaryField()[patchi][facei])/rhoAvg.boundaryField()[patchi][facei]);
            }
        }

        CpEff.correctBoundaryConditions();

    }

    volScalarField Tdiff = 0.*thermo.T();  //initialize Tdiff to 0.0

    forAll(rho.internalField(), celli)
    {
        Tdiff[celli] = Twall - thermo.T()[celli];
    }

    forAll(rho.boundaryField(), patchi)
    {
        forAll(rho.boundaryField()[patchi], facei)
        {
            Tdiff.boundaryFieldRef()[patchi][facei] = Twall - thermo.T().boundaryFieldRef()[patchi][facei];
        }
    }   
