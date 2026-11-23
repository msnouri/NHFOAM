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

Class
    Foam::RASModels::kOmegaLF

Description
    The k-omega model of Wilcox (2006) with the stabilisation of Larsen &
    Fuhrman (2018) for flows under surface waves.

        dk/dt + u.grad(k) = P - betaStar k omega
                          + div[(nu + sigmaK k/omega) grad(k)]
        domega/dt + u.grad(omega) = gamma (omega/k) P - beta omega^2
                          + (sigmaD/omega) grad(k).grad(omega)
                          + div[(nu + sigmaOmega k/omega) grad(omega)]

    with P = nut p0, p0 = 2 S:S, beta = beta0 fBeta (vortex stretching,
    Pope 1978), sigmaD = sigmaDo where grad(k).grad(omega) > 0, else 0, and

        nut = k/omegaTilde,
        omegaTilde = max[omega, lambda1 sqrt(p0/betaStar),
                         lambda2 (beta/(betaStar gamma)) (p0/pOmega) omega],

    pOmega = 2 Omega:Omega.  The lambda1 term is the stress limiter of Wilcox
    (2006); the lambda2 term (Larsen & Fuhrman 2018) keeps the model stable
    in nearly potential flow (pOmega << p0), where the standard closures
    produce turbulence under non-breaking waves.  lambda2 = 0 gives the
    standard Wilcox (2006) model.

    Default coefficients
        kOmegaLFCoeffs
        {
            betaStar    0.09;
            beta0       0.0708;
            gamma       0.52;       // 13/25
            sigmaK      0.6;        // sigma*
            sigmaOmega  0.5;        // sigma
            sigmaDo     0.125;
            lambda1     0.875;      // 7/8 (Wilcox's Clim)
            lambda2     0.05;
            vortexStretching true;  // fBeta; false: beta = beta0
        }

    References
        Wilcox, D.C. (2006). Turbulence Modeling for CFD, 3rd edn. DCW
        Industries.
        Larsen, B.E., Fuhrman, D.R. (2018). On the over-production of
        turbulence beneath surface waves in Reynolds-averaged Navier-Stokes
        models. J. Fluid Mech. 853, 419-460.

SourceFiles
    kOmegaLF.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaLF_H
#define kOmegaLF_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaLF Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaLF
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kOmegaLF(const kOmegaLF&) = delete;

        //- No copy assignment
        void operator=(const kOmegaLF&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar beta0_;
            dimensionedScalar gamma_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;
            dimensionedScalar sigmaDo_;
            dimensionedScalar lambda1_;
            dimensionedScalar lambda2_;
            Switch vortexStretching_;


        // Fields

            volScalarField k_;
            volScalarField omega_;


    // Protected Member Functions

        //- beta = beta0 fBeta
        tmp<volScalarField::Internal> beta
        (
            const volTensorField& gradU
        ) const;

        //- Limited omega of the eddy viscosity
        tmp<volScalarField> omegaTilde(const volTensorField& gradU) const;

        //- Update nut from the velocity gradient
        void correctNut(const volTensorField& gradU);

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaLF");


    // Constructors

        //- Construct from components
        kOmegaLF
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
    virtual ~kOmegaLF() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Effective diffusivity for k (Wilcox 2006: k/omega, not nut)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    sigmaK_*k_/omega_ + this->nu()
                )
            );
        }

        //- Effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    sigmaOmega_*k_/omega_ + this->nu()
                )
            );
        }

        //- Turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaLF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
