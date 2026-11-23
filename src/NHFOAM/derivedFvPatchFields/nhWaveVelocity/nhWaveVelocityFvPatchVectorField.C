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

#include "nhWaveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "gravityMeshObject.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nhWaveVelocityFvPatchVectorField::nhWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    waveDict_(),
    activeAbsorption_(true),
    absorptionMode_("shallow"),
    absorptionCelerity_(0),
    swlGiven_(false),
    swl_(0),
    etaName_("eta"),
    wavePtr_(nullptr),
    gMag_(0)
{}


Foam::nhWaveVelocityFvPatchVectorField::nhWaveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    waveDict_(dict.subDict("wave")),
    activeAbsorption_(dict.getOrDefault<bool>("activeAbsorption", true)),
    absorptionMode_("shallow"),
    absorptionCelerity_(0),
    swlGiven_(dict.found("stillWaterLevel")),
    swl_(dict.getOrDefault<scalar>("stillWaterLevel", 0)),
    etaName_(dict.getOrDefault<word>("eta", "eta")),
    wavePtr_(nullptr),
    gMag_(0)
{
    if (dict.found("absorptionCelerity"))
    {
        ITstream& is = dict.lookup("absorptionCelerity");
        const token tok(is);

        if (tok.isWord())
        {
            absorptionMode_ = tok.wordToken();
        }
        else if (tok.isNumber())
        {
            absorptionMode_ = "value";
            absorptionCelerity_ = tok.number();
        }

        if
        (
            absorptionMode_ != "shallow"
         && absorptionMode_ != "wave"
         && absorptionMode_ != "value"
        )
        {
            FatalIOErrorInFunction(dict)
                << "absorptionCelerity " << absorptionMode_ << " unknown: use"
                << " shallow, wave or a value [m/s]" << exit(FatalIOError);
        }
    }

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(vector::zero);
    }
}


Foam::nhWaveVelocityFvPatchVectorField::nhWaveVelocityFvPatchVectorField
(
    const nhWaveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    waveDict_(ptf.waveDict_),
    activeAbsorption_(ptf.activeAbsorption_),
    absorptionMode_(ptf.absorptionMode_),
    absorptionCelerity_(ptf.absorptionCelerity_),
    swlGiven_(ptf.swlGiven_),
    swl_(ptf.swl_),
    etaName_(ptf.etaName_),
    wavePtr_(nullptr),
    gMag_(0)
{}


Foam::nhWaveVelocityFvPatchVectorField::nhWaveVelocityFvPatchVectorField
(
    const nhWaveVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    waveDict_(ptf.waveDict_),
    activeAbsorption_(ptf.activeAbsorption_),
    absorptionMode_(ptf.absorptionMode_),
    absorptionCelerity_(ptf.absorptionCelerity_),
    swlGiven_(ptf.swlGiven_),
    swl_(ptf.swl_),
    etaName_(ptf.etaName_),
    wavePtr_(nullptr),
    gMag_(0)
{}


Foam::nhWaveVelocityFvPatchVectorField::nhWaveVelocityFvPatchVectorField
(
    const nhWaveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    waveDict_(ptf.waveDict_),
    activeAbsorption_(ptf.activeAbsorption_),
    absorptionMode_(ptf.absorptionMode_),
    absorptionCelerity_(ptf.absorptionCelerity_),
    swlGiven_(ptf.swlGiven_),
    swl_(ptf.swl_),
    etaName_(ptf.etaName_),
    wavePtr_(nullptr),
    gMag_(0)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::regularWave&
Foam::nhWaveVelocityFvPatchVectorField::wave() const
{
    if (!wavePtr_)
    {
        gMag_ = mag(meshObjects::gravity::New(db().time()).value());

        if (!swlGiven_ && db().foundObject<IOdictionary>("NHFOAMProperties"))
        {
            swl_ =
                db().lookupObject<IOdictionary>("NHFOAMProperties")
               .getOrDefault<scalar>("stillWaterLevel", 0);
        }

        Info<< "nhWaveVelocity on patch " << patch().name() << ":" << endl;
        wavePtr_.reset(new regularWave(waveDict_, gMag_, swl_));
    }

    return *wavePtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nhWaveVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const regularWave& w = wave();
    const scalar t = db().time().value();

    const vectorField& Cf = patch().Cf();

    vectorField Up(Cf.size());
    forAll(Cf, i)
    {
        Up[i] = w.velocity(Cf[i], t);
    }

    if (activeAbsorption_ && db().foundObject<volScalarField>(etaName_))
    {
        // Elevation of the boundary columns at the old time level against
        // the target at their centres
        const volScalarField& eta =
            db().lookupObject<volScalarField>(etaName_);
        const scalarField etaC
        (
            eta.boundaryField()[patch().index()].patchInternalField()
        );
        const vectorField Cn(patch().Cn());
        const vectorField n(patch().nf());
        const scalar t0 = t - db().time().deltaTValue();

        const scalar h = w.depth();
        scalar c = Foam::sqrt(gMag_*h);
        if (absorptionMode_ == "wave")
        {
            c = w.omega()/w.k();
        }
        else if (absorptionMode_ == "value")
        {
            c = absorptionCelerity_;
        }

        forAll(Up, i)
        {
            const scalar etaR = etaC[i] - swl_ - w.eta(Cn[i], t0);
            Up[i] += (c/h)*etaR*n[i];
        }
    }

    vectorField::operator=(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::nhWaveVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("activeAbsorption", activeAbsorption_);
    if (absorptionMode_ == "value")
    {
        os.writeEntry("absorptionCelerity", absorptionCelerity_);
    }
    else
    {
        os.writeEntry("absorptionCelerity", absorptionMode_);
    }
    if (swlGiven_)
    {
        os.writeEntry("stillWaterLevel", swl_);
    }
    os.writeEntryIfDifferent<word>("eta", "eta", etaName_);
    waveDict_.writeEntry("wave", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        nhWaveVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
