{
    volScalarField& he = thermo.he();

    if (useOverallU || porousThermalConductivity)
    {
        #include "preCalcEE.H"
    }


    fvScalarMatrix EEqn  // warning: this approach is intended for steady state / LTS use only
    (
        porosity*fvm::ddt(rho, he)
      + (
            mvConvection.valid()
          ? mvConvection->fvmDiv(phi, he)
          : fvm::div(phi, he)
        )
      + porosity*fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection.valid()
            ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
            : fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      - fvm::laplacian(alphaEff, he)
//      + thermophysicalTransport.divq(he)
     ==
        rho*(U&g)
      + porosity*reaction.Qdot()
      + porosity*fvModels.source(rho, he)
    );

    if (useOverallU)
    {
          Info << "Overall U (heat transfer coefficient) used in EEqn" << endl;
          EEqn += - overallU*Tdiff*4/(diamTube*UnitLength);  // U W/m2-K;  Tdiff K;   dH m
    }

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
