/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::RASModels::multiphaseMixtureERzetaF

Group
    grpRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible and compressible
    flows including rapid distortion theory (RDT) based compression term.

    Reference:
    \verbatim
        Standard model:
            Launder, B. E., & Spalding, D. B. (1972).
            Lectures in mathematical models of turbulence.

            Launder, B. E., & Spalding, D. B. (1974).
            The numerical computation of turbulent flows.
            Computer methods in applied mechanics and engineering,
            3(2), 269-289.

        For the RDT-based compression term:
            El Tahry, S. H. (1983).
            k-epsilon equation for compressible reciprocating engine flows.
            Journal of Energy, 7(4), 345-353.
    \endverbatim

    The default model coefficients are
    \verbatim
        multiphaseMixtureERzetaFCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    multiphaseMixtureERzetaF.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMixtureERzetaF_H
#define multiphaseMixtureERzetaF_H

#include "RASModel.H"
#include "multiphaseMixtureEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class multiphaseMixtureERzetaF Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class multiphaseMixtureERzetaF
:
    public multiphaseMixtureEddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

    	mutable multiphaseMixtureERzetaF<BasicTurbulenceModel> *firstPhaseTurbulencePtr_;

        //- No copy construct
        multiphaseMixtureERzetaF(const multiphaseMixtureERzetaF&) = delete;

        //- No copy assignment
        void operator=(const multiphaseMixtureERzetaF&) = delete;

        //- Return the turbulence model for the other phase
        multiphaseMixtureERzetaF<BasicTurbulenceModel> & firstPhaseTurbulence() const;

        bool initialized;


protected:

    // Protected data

        // Model coefficients
        dimensionedScalar Cmu_;
        dimensionedScalar COmega2_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar sigmaK_;
        dimensionedScalar sigmaOmega_;
        dimensionedScalar sigmaCDv_;
        dimensionedScalar sigmaCDt_;
        dimensionedScalar sigmaZeta_;
        dimensionedScalar CTau_;
        dimensionedScalar CL_;
        dimensionedScalar CEta_;
        dimensionedScalar Csas_;
        dimensionedScalar CT2_;
        dimensionedScalar Clim_;

      	// Fields
       	volScalarField k_;
       	volScalarField omega_;
       	volScalarField epsilon_;
       	volScalarField zeta_;
       	volScalarField f_;
        volScalarField vt_;
       	const volScalarField& yr_;
        PtrList<multiphaseMixtureERzetaF<BasicTurbulenceModel>> turbulencePtrList_;

        //- Lower limit of zeta
        dimensionedScalar mTSmall_;
        dimensionedScalar zetaMin_;
        //- Lower limit of f
        dimensionedScalar fMin_;
       
    // Protected Member Functions

        void initilizePtrList();

        //- Return time scale, Ts
        tmp<volScalarField> Tau() const;
        //- Return length scale, L
        tmp<volScalarField> L() const;

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("multiphaseMixtureERzetaF");


    // Constructors

        //- Construct from components
        multiphaseMixtureERzetaF
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~multiphaseMixtureERzetaF() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

            //- Return the effective diffusivity for omega
            tmp<volScalarField> DomegaEff() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DomegaEff",
                        this->nut_/sigmaOmega_ + this->nu()
                    )
                );
            }

            //
            tmp<volScalarField> DzetaEff() const
            {
                return tmp<volScalarField>
                (
                    new volScalarField
                    (
                        "DzetaEff",
                        this->nut_/sigmaZeta_ + this->nu()
                    )
                );
            }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        tmp<volScalarField> kk() 
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiphaseMixtureERzetaF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
