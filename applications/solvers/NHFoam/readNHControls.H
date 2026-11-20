// Move the mesh in every outer corrector (true) or only once per time step,
// to the extrapolated free surface (false, default)
const bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault<bool>("moveMeshOuterCorrectors", false)
);

// Rhie-Chow type time-derivative flux correction (as pimpleFoam's ddtCorr)
const bool ddtCorr
(
    pimple.dict().getOrDefault<bool>("ddtCorr", true)
);
