/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | NHFOAM: non-hydrostatic free-surface wave solver
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of NHFOAM, built on OpenFOAM.

    NHFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "kOmegaLF.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> kOmegaLF<BasicTurbulenceModel>::beta
(
    const volTensorField& gradU
) const
{
    tmp<volScalarField::Internal> tbeta
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::groupName("kOmegaLF:beta", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            beta0_
        )
    );

    if (vortexStretching_)
    {
        // fBeta = (1 + 85 chi)/(1 + 100 chi),
        // chi = |Omega_ij Omega_jk S_ki|/(betaStar omega)^3
        const volSymmTensorField::Internal S(symm(gradU.v()));
        const volTensorField::Internal W(skew(gradU.v()));

        const volScalarField::Internal chi
        (
            mag((W & W) && S)/pow3(betaStar_*omega_.v())
        );

        tbeta.ref() = beta0_*(1 + 85*chi)/(1 + 100*chi);
    }

    return tbeta;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLF<BasicTurbulenceModel>::omegaTilde
(
    const volTensorField& gradU
) const
{
    const volScalarField p0(2*magSqr(symm(gradU)));

    tmp<volScalarField> tw
    (
        new volScalarField
        (
            IOobject::groupName("omegaTilde", this->alphaRhoPhi_.group()),
            max(omega_, lambda1_*sqrt(p0/betaStar_))
        )
    );

    if (lambda2_.value() > 0)
    {
        const volScalarField pOmega(2*magSqr(skew(gradU)));

        volScalarField betaF
        (
            IOobject
            (
                "kOmegaLF:betaF",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            beta0_,
            calculatedFvPatchScalarField::typeName
        );
        betaF.primitiveFieldRef() = beta(gradU)().field();
        betaF.correctBoundaryConditions();

        tw.ref() = max
        (
            tw(),
            lambda2_*(betaF/(betaStar_*gamma_))
           *(p0/max(pOmega, dimensionedScalar("small", pOmega.dimensions(), SMALL)))
           *omega_
        );
    }

    return tw;
}


template<class BasicTurbulenceModel>
void kOmegaLF<BasicTurbulenceModel>::correctNut(const volTensorField& gradU)
{
    this->nut_ = k_/omegaTilde(gradU);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaLF<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaLF<BasicTurbulenceModel>::kOmegaLF
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar", this->coeffDict_, 0.09
        )
    ),
    beta0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta0", this->coeffDict_, 0.0708
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma", this->coeffDict_, 0.52
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK", this->coeffDict_, 0.6
        )
    ),
    sigmaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaOmega", this->coeffDict_, 0.5
        )
    ),
    sigmaDo_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaDo", this->coeffDict_, 0.125
        )
    ),
    lambda1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambda1", this->coeffDict_, 0.875
        )
    ),
    lambda2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambda2", this->coeffDict_, 0.05
        )
    ),
    vortexStretching_
    (
        Switch::lookupOrAddToDict
        (
            "vortexStretching", this->coeffDict_, true
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaLF<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaDo_.readIfPresent(this->coeffDict());
        lambda1_.readIfPresent(this->coeffDict());
        lambda2_.readIfPresent(this->coeffDict());
        vortexStretching_.readIfPresent("vortexStretching", this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmegaLF<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();

    const volScalarField::Internal GbyNu
    (
        gradU.v() && dev(twoSymm(gradU.v()))
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);

    // Production of omega, gamma (omega/k) P = gamma (omega/omegaTilde) p0
    const volScalarField::Internal omegaRatio
    (
        omega_()/omegaTilde(gradU)()()
    );

    const volScalarField::Internal betaI(beta(gradU));

    // Cross diffusion (sigmaD/omega) grad(k).grad(omega), where positive
    const volScalarField::Internal CDkOmega
    (
        sigmaDo_*max
        (
            (fvc::grad(k_) & fvc::grad(omega_))().v(),
            dimensionedScalar
            (
                "zero",
                k_.dimensions()*omega_.dimensions()/sqr(dimLength),
                0
            )
        )/omega_()
    );

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha()*rho()*GbyNu*omegaRatio
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(betaI*alpha()*rho()*omega_(), omega_)
      + alpha()*rho()*CDkOmega
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    // Turbulence kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(betaStar_*alpha()*rho()*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
