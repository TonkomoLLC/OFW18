/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::heRhoThermoPorous

Description
    Energy for a mixture based on density. Optionally calculates effective thermal 
    conducitivty and effective heat capacity for a speciifed porous media region.

    Based on the references:
        Zehner P. and Schlünder E.U, "Thermal Conductivity of Granular Materials 
        at Moderate Temperatures", Chemie.Ingr. Tech., 42, 1970, 933-941.

        Specchia V., Baldi G., Sicardi S., "Heat transfer in-packed bed reactors 
        with one phase flow", Chem. Engng Commun, 4, 1980, 361-380.
        
SourceFiles
    heRhoThermoPorous.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoThermoPorous_H
#define heRhoThermoPorous_H

#include "rhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoThermoPorous Declaration
\*---------------------------------------------------------------------------*/

template<class BasicRhoThermo, class MixtureType>
class heRhoThermoPorous
:
    public heThermo<BasicRhoThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


public:

    //- Runtime type information
    TypeName("heRhoThermoPorous");


    // Constructors

        //- Construct from mesh and phase name
        heRhoThermoPorous
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        heRhoThermoPorous(const heRhoThermoPorous<BasicRhoThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~heRhoThermoPorous();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const heRhoThermoPorous<BasicRhoThermo, MixtureType>&
        ) = delete;

    // - Porous media thermal conductivity

        // - Switch to use porous media thermal conductivity 
        Switch usePorousKappa;

        // - List of cells in the porous cell zone
        labelList cellsZone;

        // - Average porosity of porous cell zone
        scalar porosity;

        // - Solid thermal conductivity, W/m/K
        scalar kappaSolid;

        // - Solid heat capacity, J/kg/K
        scalar cPSolid;

        // - Solid density, kg/m3
        scalar rhoSolid;

        // - Solid particle diameter, m
        scalar diamParticle;

        // - Superficial velocity, m/s
        scalar USuperficial;

        // - Tube diameter, m
        scalar diamTube;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heRhoThermoPorous.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
