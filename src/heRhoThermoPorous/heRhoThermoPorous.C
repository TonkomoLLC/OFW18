/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "heRhoThermoPorous.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::heRhoThermoPorous<BasicRhoThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli]);

        //const fvMesh& mesh_ = this->T_.mesh();       
        //const volVectorField& USuperficial = mesh_.lookupObject<volVectorField>("U");  //hard coded for velocity field name! 

        if (usePorousKappa)
        {  
            forAll(cellsZone, cI)
            {
                if (cellsZone[cI] == celli)
                {
                    // Effective radial solid thermal conductivity, Zhener and Schlunder 
                    scalar kfks =  kappaCells[celli] / kappaSolid;
                    scalar kskf = kappaSolid / kappaCells[celli];
                    scalar Cterm = 1.25;      // setup for spheres!
                    scalar Bterminner = (1-porosity)/porosity;
                    scalar Bterm = Cterm*pow(Bterminner,10./9.);
                    scalar Bterm2 = 1 - kfks*Bterm;
                    scalar Bterm3 = (1 - kfks)*Bterm;
                    scalar term0 = 1 - sqrt(1-porosity); 
                    scalar term1 = sqrt(1-porosity);
                    scalar term2 = 2.0/Bterm2;
                    scalar term3 = Bterm3 / pow(Bterm2,2.0);
                    scalar term4 = log(kskf/Bterm);
                    scalar term5 = (Bterm+1)/2;
                    scalar term6 = (Bterm - 1) / Bterm2;
                    scalar kappaEffSolid = kappaCells[celli]*(term0 + term1*(term2*(term3*term4-term5-term6)));
                       
                    // Effective radial fluid thermal conductivity, Specchia et al.
                    //scalar ReNum = thermoMixture.rho(pCells[celli], TCells[celli])*mag(USuperficial[celli])*diamParticle/muCells[celli];
                    scalar ReNum = thermoMixture.rho(pCells[celli], TCells[celli])*USuperficial*diamParticle/muCells[celli];
                    scalar PrNum = CpCells[celli] * muCells[celli]/kappaCells[celli];
                    scalar Re_a = ReNum; 
                    scalar Nnum = diamParticle/diamTube;
                    scalar B_Specchia = 1 + 19.4 * pow(Nnum,2.0);
                    scalar Pe_H = 8.65*B_Specchia;
                    scalar kappaEffConvection = kappaCells[celli]*(Re_a*PrNum/Pe_H);

                    kappaCells[celli] = kappaEffSolid + kappaEffConvection;

                    // Effective heat capacity
                    scalar rhoCelli = thermoMixture.rho(pCells[celli], TCells[celli]);
                    scalar rhoAverageCelli = porosity*rhoCelli + (1 - porosity)*rhoSolid;
                    scalar cPCelli = ((CpCells[celli] * porosity * rhoCelli) 
                                + ((1 - porosity)*cPSolid*rhoSolid))/rhoAverageCelli;
                    scalar cPFactor = CpCells[celli] / cPCelli; 

                    kappaCells[celli] *= cPFactor; 
                }
            }
        }
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::heRhoThermoPorous<BasicRhoThermo, MixtureType>::heRhoThermoPorous
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicRhoThermo, MixtureType>(mesh, phaseName)
{

    // Setup porous properties
    IOdictionary porousPropertiesDict
    (
        IOobject
        (
            "porousProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Read properties for effective thermal conductivity calculation
    usePorousKappa = porousPropertiesDict.lookupOrDefault<Switch>("porousThermalConductivity", false);
    diamParticle = porousPropertiesDict.lookupOrDefault<scalar>("diamParticle",0.005);    
    USuperficial = porousPropertiesDict.lookupOrDefault<scalar>("USuperficial",0.005);   // not used; get U from field data
    diamTube = porousPropertiesDict.lookupOrDefault<scalar>("diamTube",0.05);    
    porosity = porousPropertiesDict.lookupOrDefault<scalar>("porosity",1);
    kappaSolid = porousPropertiesDict.lookupOrDefault<scalar>("kappaSolid",20);
    cPSolid = porousPropertiesDict.lookupOrDefault<scalar>("cPSolid",900);
    rhoSolid = porousPropertiesDict.lookupOrDefault<scalar>("rhoSolid",2000);
    word porousZoneName = porousPropertiesDict.lookupOrDefault<word>("porousZoneName","porosity");

    label zoneID = mesh.cellZones().findZoneID(porousZoneName);

    if ((zoneID != -1) && usePorousKappa)
    {
        Info << "Porosity zone ID for thermal conductivity adjustments = " << zoneID << endl;
        Info << "Porosity, kappaSolid, porousZoneName : " << porosity << " " << kappaSolid << " " << porousZoneName << endl;

        const cellZone& cellsPorosity = mesh.cellZones()[zoneID];
        const labelList& cellLabels = mesh.cellZones()[zoneID];
        const meshCellZones& zoneMesh = cellsPorosity.meshZones();
        cellsZone = zoneMesh[zoneID];   //list of all porosityZone cell ID's
        label porositySize = returnReduce(cellLabels.size(), sumOp<label>());
        Info << "Number of cells in porous zone: ";
        Info << porositySize << endl;
    }

    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::heRhoThermoPorous<BasicRhoThermo, MixtureType>::~heRhoThermoPorous()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::heRhoThermoPorous<BasicRhoThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
