/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "multiphaseMixtureERzetaF.H"
#include "multiphaseSystem.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallFvPatch.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> multiphaseMixtureERzetaF<BasicTurbulenceModel>::Tau() const
{
    return 1.0/omega_;
}
template<class BasicTurbulenceModel>
tmp<volScalarField> multiphaseMixtureERzetaF<BasicTurbulenceModel>::L() const
{
    return CL_*pow(k_,0.5)/omega_;
}

template<class BasicTurbulenceModel>
void multiphaseMixtureERzetaF<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = this->Cmu_*zeta_*k_*Tau();
    this->vt_ = this->nut_;
}

template<class BasicTurbulenceModel>
void multiphaseMixtureERzetaF<BasicTurbulenceModel>::initilizePtrList()
{

  Info << "--------------------------------" << endl;
  Info << "Initializing phaseModelsPtrList:" << endl;
  const volVectorField& U = this->U_;
  const transportModel& currPhase = this->transport();
  const multiphaseSystem& fluid = refCast<const multiphaseSystem>(currPhase.fluid());

  forAll(fluid.phases(),phasei)
  {
    Info <<"    Appending: " << fluid.phases()[phasei].name() << endl;

    turbulencePtrList_.append(

        &const_cast<multiphaseMixtureERzetaF<BasicTurbulenceModel>&>
        (
                U.db().lookupObject<multiphaseMixtureERzetaF<BasicTurbulenceModel>>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        fluid.phases()[phasei].name()
                    )
                )
        )
    );
  }
  Info << "--------------------------------" << endl;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
multiphaseMixtureERzetaF<BasicTurbulenceModel>::multiphaseMixtureERzetaF
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    multiphaseMixtureEddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    firstPhaseTurbulencePtr_(nullptr),

    initialized(false),

    Cmu_    (dimensioned<scalar>::lookupOrAddToDict("Cmu",          this->coeffDict_, 0.22)),
    COmega2_        (dimensioned<scalar>::lookupOrAddToDict("COmega2",      this->coeffDict_, 0.9)),
    C1_         (dimensioned<scalar>::lookupOrAddToDict("C1",           this->coeffDict_, 0.4)),
    C2_         (dimensioned<scalar>::lookupOrAddToDict("C2",           this->coeffDict_, 0.65)),
    sigmaK_     (dimensioned<scalar>::lookupOrAddToDict("sigmaK",       this->coeffDict_, 1.1)),
    sigmaOmega_     (dimensioned<scalar>::lookupOrAddToDict("sigmaOmega",       this->coeffDict_, 1.1)),
    sigmaCDv_       (dimensioned<scalar>::lookupOrAddToDict("sigmaCDv",         this->coeffDict_, 1.2)),
    sigmaCDt_       (dimensioned<scalar>::lookupOrAddToDict("sigmaCDt",         this->coeffDict_, 1.6)),
    sigmaZeta_      (dimensioned<scalar>::lookupOrAddToDict("sigmaZeta",        this->coeffDict_, 1.2)),
    CTau_       (dimensioned<scalar>::lookupOrAddToDict("CTau",         this->coeffDict_, 6.0)),
    CL_         (dimensioned<scalar>::lookupOrAddToDict("CL",           this->coeffDict_, 0.36)),
    CEta_       (dimensioned<scalar>::lookupOrAddToDict("CEta",         this->coeffDict_, 85)),
    Csas_       (dimensioned<scalar>::lookupOrAddToDict("Csas",         this->coeffDict_, 4)),
    CT2_        (dimensioned<scalar>::lookupOrAddToDict("CT2",          this->coeffDict_, 1.0)),
    Clim_       (dimensioned<scalar>::lookupOrAddToDict("Clim",         this->coeffDict_, 0.0)),

    k_      (IOobject ("k"      , this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE ), this->mesh_),
    omega_  (IOobject ("omega"  , this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE ), this->mesh_),
    epsilon_(IOobject ("epsilon", this->runTime_.timeName(),this->mesh_, IOobject::NO_READ,   IOobject::NO_WRITE ), omega_*k_  ),
    zeta_   (IOobject ("zeta"   , this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE ), this->mesh_),
    f_      (IOobject ("f"      , this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE ), this->mesh_),
    vt_     (IOobject ("vt"     , this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE ), this->mesh_),
    yr_(wallDist::New(this->mesh_).y()),
    turbulencePtrList_ (),
    mTSmall_("mTSmall", dimensionSet(0, 0, -1, 0, 0, 0, 0),1e-10),
    zetaMin_("zetaMin", dimless, 1e-10),
    fMin_("fMin", dimless/dimTime, 1e-10)

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool multiphaseMixtureERzetaF<BasicTurbulenceModel>::read()
{
    if (multiphaseMixtureEddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        COmega2_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaCDv_.readIfPresent(this->coeffDict());
        sigmaCDt_.readIfPresent(this->coeffDict());
        sigmaZeta_.readIfPresent(this->coeffDict());
        CTau_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        CEta_.readIfPresent(this->coeffDict());
        Csas_.readIfPresent(this->coeffDict());
        CT2_.readIfPresent(this->coeffDict());
        Clim_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
multiphaseMixtureERzetaF<BasicTurbulenceModel>&
multiphaseMixtureERzetaF<BasicTurbulenceModel>::firstPhaseTurbulence() const
{
    if (!firstPhaseTurbulencePtr_) // If the pointer points to NULL
    {
        const volVectorField& U = this->U_;
        const transportModel& currPhase = this->transport();
        const multiphaseSystem& fluid = refCast<const multiphaseSystem>(currPhase.fluid());

        const transportModel& firstPhase = fluid.phases()[0];

        firstPhaseTurbulencePtr_ =

        &const_cast<multiphaseMixtureERzetaF<BasicTurbulenceModel>&>
            (
                U.db().lookupObject<multiphaseMixtureERzetaF<BasicTurbulenceModel>>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        firstPhase.name()
                    )
                )
            );
    }

    return *firstPhaseTurbulencePtr_;
}



template<class BasicTurbulenceModel>
void multiphaseMixtureERzetaF<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    const transportModel& currPhase = this->transport();
    const multiphaseSystem& fluid = refCast<const multiphaseSystem>(currPhase.fluid());
    
    this->nut_ = this->vt_;

    // Only solve the mixture turbulence for the first phase
    if (&currPhase != &fluid.phases()[0])
    {
        //Info << "Assigning turbulence fields for phase: " << currPhase.name() << endl;
        Info << "Phase: " << currPhase.name() << " is not the first phase, skipping" << endl;
        this->firstPhaseTurbulence();
        return;
    }
    Info << "Solving turbulence for first phase: " << currPhase.name() << endl;

    if(!initialized)
    {
        initilizePtrList();
        initialized = true;
    }
//***************************************************************************
// Local references 

   const surfaceScalarField& phim = fluid.phi(); // mixture volumetric flux
   const volVectorField U(fluid.U());

   volScalarField& nut = this->nut_;
   fv::options& fvOptions(fv::options::New(this->mesh_));
   volScalarField divU(fvc::div(fvc::absolute(phim, U)));
   const volTensorField gradU(fvc::grad(U));

   const volScalarField S2(2*magSqr(dev(symm(gradU))));
   const volScalarField G(this->GName(), nut*S2);
   const volScalarField T(this->Tau());
   const volScalarField L(this->L());

//******************************************************************************************
   const volScalarField COmega1
   (
     "COmega1",
     (1.4*(1.0 + (0.012/(zeta_+zetaMin_))))-1
   );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAS TERM LIMITER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    dimensionedScalar Psaslim("Psaslim", dimensionSet(0, 0, -2, 0, 0), 0.0);
    volScalarField T1_ = 40.0*1.775*0.41*mag(fvc::laplacian(U))*sqrt(k_);
    volScalarField T2_ = 3.0*k_*max(pow(omega_,-2.0)*(fvc::grad(omega_) & fvc::grad(omega_)), pow(k_,-2.0)*(fvc::grad(k_) & fvc::grad(k_)));
    volScalarField Psas = max(0.003*(T1_ - CT2_*T2_), Psaslim);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    volScalarField CDv("CDv", (2.0/k_*this->nu()/sigmaCDv_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CDt("CDt", (2.0/k_*nut/sigmaCDt_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CD("CD", CDv+max(CDt,dimensionedScalar("0.0", dimless/dimTime, 0.0)));


// *********************************************************************************************************************************************************************//
// Solve for turbulence fields

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phim, omega_)
      + fvm::SuSp(-fvc::div(phim), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
      ==
        COmega1*G/k_*omega_
      - fvm::SuSp(COmega2_*omega_, omega_)
      + fvm::Sp(CD, omega_)
      + Csas_*Psas
    );  

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());

    #include "wallTreatmentOmega.H"

    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phim, k_)
      + fvm::SuSp(-fvc::div(phim), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(omega_, k_)
      + fvOptions(k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L), f_)
      - (C1_+(C2_*G/(omega_*k_)))*((zeta_ - 2.0/3.0))/(sqr(L)*T)
    );

    fEqn.ref().relax();
    fvOptions.constrain(fEqn.ref());
    solve(fEqn);
    fvOptions.correct(f_);
    bound(f_, fMin_);

    volScalarField fTilda = f_ - 2.0*this->nu()*zeta_/sqr(yr_);

    // zeta Equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(zeta_)
      + fvm::div(phim, zeta_)
      + fvm::SuSp(-fvc::div(phim), zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
      ==
        fTilda
      - fvm::Sp(G/k_, zeta_)
      + fvOptions(zeta_)
    );

    zetaEqn.ref().relax();
    fvOptions.constrain(zetaEqn.ref());
    solve(zetaEqn);
    fvOptions.correct(zeta_);
    bound(zeta_, zetaMin_);
    zeta_ = min(zeta_, 2.0);

    epsilon_ = omega_*k_;


    correctNut();

// *********************************************************************************************************************************************************************//
// Assign turbulence fields and write

    
    
    for(int i = 1;i< turbulencePtrList_.size();i++)
    {
      Info << "Assigning turbulence fields to phase " << fluid.phases()[i].name() <<" :";
      turbulencePtrList_[i].k_       = k_;
      Info <<" k,";
      turbulencePtrList_[i].omega_   = omega_;
      Info <<" omega,";
      turbulencePtrList_[i].epsilon_ = epsilon_;
      Info <<" epsilon,";
      turbulencePtrList_[i].zeta_   = zeta_;
      Info <<" zeta,";
      turbulencePtrList_[i].f_       = f_;
      Info <<" f,";
      turbulencePtrList_[i].vt_      = this->vt_;
      Info <<" vt"<<endl;

    }

    if(this->runTime_.write())
    {
        k_.write();
        omega_.write();
        epsilon_.write();
        zeta_.write();
        f_.write();
        vt_.write();
    }

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
