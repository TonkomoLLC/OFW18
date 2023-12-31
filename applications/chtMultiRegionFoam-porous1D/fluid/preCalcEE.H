    alphaEff = thermo.kappa()/thermo.Cp();
    
    if(useOverallU || porousThermalConductivity)
    {
        // Calculate transfer numbers
        volScalarField Repsi = (rho * (mag(U)+small*UnitMeterPerSecond) * (diamParticle * UnitLength))/ (thermo.mu() * porosity); 
        volScalarField Rep = (rho * (mag(U)+small*UnitMeterPerSecond) * (diamParticle * UnitLength))/ thermo.mu(); 
        volScalarField Ret = (rho * (mag(U)+small*UnitMeterPerSecond) * (diamTube * UnitLength))/ thermo.mu(); 
        volScalarField PrNum = (thermo.Cp() * thermo.mu()) / thermo.kappa(); 

        if (useUSuperficial)
        {
            Repsi = (rho * USuperficial * UnitMeterPerSecond * (diamParticle * UnitLength))/ (thermo.mu() * porosity);  
            Rep = (rho * USuperficial * UnitMeterPerSecond * (diamParticle * UnitLength))/ thermo.mu();  
            Ret = (rho * USuperficial * UnitMeterPerSecond * (diamTube * UnitLength))/ thermo.mu();
        }

        scalar Nfactor = diamTube/diamParticle;  

        // effective radial solid conductivity, P. Zehner and E.U. Schlünder Chemie. Ingr. Tech., 42 (1970), 933.
        forAll(rho.internalField(), celli)
        {
            scalar kfks = thermo.kappa()[celli] / kappaSolid;
            scalar kskf = kappaSolid / thermo.kappa()[celli];
            scalar Cterm = 1.25;      // setup for spheres!
            scalar Bterminner = (1.-porosity+small)/porosity;
            scalar Bterm = Cterm*std::pow(Bterminner,10./9.);
            scalar Bterm2 = 1. - kfks*Bterm;
            scalar Bterm3 = (1. - kfks)*Bterm;
            scalar term0 = 1. - std::sqrt(1.-porosity+small); 
            scalar term1 = std::sqrt(1.-porosity+small);
            scalar term2 = 2.0/Bterm2;
            scalar term3 = Bterm3 / std::pow(Bterm2,2.0);
            scalar term4 = std::log(kskf/Bterm);
            scalar term5 = (Bterm+1.)/2.;
            scalar term6 = (Bterm - 1.) / Bterm2;

            scalar Nnum = diamParticle/diamTube;
            scalar B_Specchia = 1. + 19.4 * std::pow(Nnum,2.0);
            scalar Pe_H = 8.65*B_Specchia;
            scalar kappaEffConvection = thermo.kappa()[celli]*(Rep[celli]*PrNum[celli]/Pe_H);

            kappaEffSolid[celli] = thermo.kappa()[celli]*(term0 + term1*(term2*(term3*term4-term5-term6)));
            kappaEff[celli] = kappaEffConvection + kappaEffSolid[celli];
        }    

        forAll(rho.boundaryField(), patchi)
        {
            forAll(rho.boundaryField()[patchi], facei)
            {
                scalar kfks = thermo.kappa().boundaryField()[patchi][facei]/ kappaSolid;
                scalar kskf = kappaSolid / thermo.kappa().boundaryField()[patchi][facei];
                scalar Cterm = 1.25;      // setup for spheres!
                scalar Bterminner = (1.-porosity+small)/porosity;
                scalar Bterm = Cterm*std::pow(Bterminner,10./9.);
                scalar Bterm2 = 1. - kfks*Bterm;
                scalar Bterm3 = (1. - kfks)*Bterm;
                scalar term0 = 1. - std::sqrt(1.-porosity+small); 
                scalar term1 = std::sqrt(1.-porosity+small);
                scalar term2 = 2.0/Bterm2;
                scalar term3 = Bterm3 / std::pow(Bterm2,2.0);
                scalar term4 = std::log(kskf/Bterm);
                scalar term5 = (Bterm+1.)/2.;
                scalar term6 = (Bterm - 1.) / Bterm2;

                scalar Nnum = diamParticle/diamTube;
                scalar B_Specchia = 1. + 19.4 * std::pow(Nnum,2.0);
                scalar Pe_H = 8.65*B_Specchia;
                scalar kappaEffConvection = thermo.kappa().boundaryField()[patchi][facei]
                    *(Rep.boundaryField()[patchi][facei]*PrNum.boundaryField()[patchi][facei]/Pe_H);
                   
                kappaEffSolid.boundaryFieldRef()[patchi][facei] = 
                    thermo.kappa().boundaryField()[patchi][facei]*(term0 + term1*(term2*(term3*term4-term5-term6)));
                kappaEff.boundaryFieldRef()[patchi][facei] = kappaEffConvection + kappaEffSolid.boundaryField()[patchi][facei];
            }
        }    

        kappaEffSolid.correctBoundaryConditions();
        kappaEff.correctBoundaryConditions();

        // Calculate wall heat transfer coefficient
        // H. Martin and M. Nilles. Radiale Warmeleitung in durchstromten Schuttungsrohren. Chem. Ing. Tech. 65 (1993), 14681477.
        volScalarField Nuw = (1.3 + 5./Nfactor)*(kappaEffSolid/thermo.kappa()) + 0.19*pow(PrNum,0.33333)*pow(Rep,0.75);
        volScalarField hw = Nuw * thermo.kappa() / (diamParticle*UnitLength);


        // Bi = hw * (d_t/2) / Keff
        // A.G. Dixon. An improved equation for the overall heat transfer coefficient in packed beds. 
        // Chemical Engineering and Processing 35 (1996), 323-331.
        
        Bi = hw * (diamTube*UnitLength/2) / kappaEff;

        forAll(rho, celli)
        {
            forAll(cellsZone, cI)
            {
                if (cellsZone[cI] == celli)
                {    
                    
                    OneOverU[celli] = 1./hw[celli] + (diamTube)/(6.*kappaEff[celli]) *((Bi[celli]+3.)/(Bi[celli]+4.));  
                    if (useOverallUconst)
                        overallU[celli] = overallUconst;
                    else
                        overallU[celli] = 1. / OneOverU[celli];

                    alphaEff[celli] = kappaEff[celli] / CpEff[celli];
                }
            }
        }
        overallU.correctBoundaryConditions();
    }

    alphaEff.correctBoundaryConditions();
