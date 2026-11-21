// Move the mesh in every outer corrector (true) or only once per time step,
// to the extrapolated free surface (false, default)
const bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault<bool>("moveMeshOuterCorrectors", false)
);

// Convecting flux of the first outer corrector: extrapolated in time from
// the face velocities of levels n and n-1 (true, default with backward),
// which keeps the convection term second order in time with one outer
// corrector, or the flux of level n (false: first order in time unless
// nOuterCorrectors > 1)
const bool fluxExtrapolation
(
    pimple.dict().getOrDefault<bool>("fluxExtrapolation", true)
);

// Rhie-Chow type time-derivative flux correction (as pimpleFoam's ddtCorr)
const bool ddtCorr
(
    pimple.dict().getOrDefault<bool>("ddtCorr", true)
);
