#!/bin/bash
# Regenerate the figures of doc/figures (README and doc/report) from the
# validation runs.  The run directories below are those of the author's
# machine: set RUNS, CMP, MA, LAB and NHW to your own copies.
#
#   RUNS : NHFOAM runs (interface layout: $RUNS_IF, cell layout: $RUNS_CELL)
#   CMP  : NHWAVE and REEF3D::NHFLOW runs (other_models.py)
#   MA   : runs of the Ma et al. (2012) tests (ma2012.py)
#   LAB  : fengyanshi/BENCHMARK_FUNWAVE
#   NHW  : JimKirby/NHWAVE (lab records of the landslide benchmark)
set -e
cd "$(dirname "$0")"
RUNS_IF=${RUNS_IF:-/home/claude/runs_if}
RUNS_CELL=${RUNS_CELL:-/home/claude/runs}
CMP=${CMP:-/home/claude/runs_cmp}
MA=${MA:-/home/claude/runs_ma}
LAB=${LAB:-/home/claude/ext/BENCHMARK_FUNWAVE}
NHW=${NHW:-/home/claude/ext/NHWAVE}
OUT=${OUT:-../doc/figures}
mkdir -p "$OUT"

# Basic checks
python3 plot_basin.py $RUNS_IF/ibasin $OUT/circular_basin_3d.png
python3 plot_solitary.py $RUNS_IF/isol $OUT/solitary_wave.png

# Dispersion and layers
python3 plot_dispersion.py $RUNS_IF/disp2/dispersion.csv $OUT/dispersion_vs_layers.png \
    "Linear dispersion (interface pressure): standing waves, 64 cells per wavelength"
python3 plot_grading.py $RUNS_IF/grading/grading.csv $OUT/layer_grading.png \
    "Standing waves, interface pressure: uniform vs surface-refined \$\\sigma\$-layers"
python3 plot_dispersion_models.py $RUNS_IF/disp2/dispersion.csv $CMP/disp_nw \
    $OUT/dispersion_models.png $RUNS_CELL/dispersion.csv

# Convection
python3 plot_flume.py 1.01 $OUT/flume_convection.png \
    $RUNS_IF/flume/iflume_f20:"Gauss linearUpwind grad(U)" \
    $RUNS_IF/flume/iflume_lust:"Gauss LUST grad(U)" \
    $RUNS_IF/flume/iflume_llv:"Gauss limitedLinearV 1" \
    $RUNS_IF/flume/iflume_cdlin:"Gauss linear (central)"

# Submerged bar, 10 layers, all models over the same window
# (NHWAVE and REEF3D::NHFLOW: 50 s runs in case A, 65 s runs in case C)
for c in A C; do
    opt=""; other=""
    [ $c = C ] && { opt="--caseC"; other=65; }
    python3 plot_bar.py $opt --tend 50 $LAB/car_luth_AC/lab_data/Case$c $OUT/bar_models_$c.png \
        $RUNS_IF/ibar${c}_L10:NHFOAM \
        $RUNS_CELL/bar${c}_l65:"NHFOAM, cell-centred pressure" \
        nhwave:$CMP/nw_bar${c}$other:NHWAVE \
        reef3d:$CMP/r3_bar${c}$other:"REEF3D::NHFLOW"
done

# Seabed motion
python3 plot_uplift.py $RUNS_IF/iuplift $OUT/bottom_uplift.png
python3 plot_enet_grilli.py $NHW/examples/landslide_benchmarks/BM2_enet_grilli/d61 \
    $OUT/landslide_models.png \
    $RUNS_IF/ieg_llv2:NHFOAM \
    $RUNS_CELL/eg_llv:"NHFOAM, cell-centred pressure" \
    nhwave:$CMP/nw_eg:NHWAVE

# Breaking
python3 plot_hansen_svendsen.py $LAB/car_hansen_svendsen/postprocessing/ht_setup_061071.txt 1.667 \
    $OUT/hansen_svendsen_models.png \
    $RUNS_IF/ihs061_or:NHFOAM \
    $RUNS_CELL/hs061:"NHFOAM, cell-centred pressure" \
    reef3d:$CMP/r3_hs:"REEF3D::NHFLOW"

# Ma et al. (2012) tests against NHWAVE
python3 plot_ma2012_basin.py $MA/basin $OUT/ma2012_basin.png
python3 plot_ma2012_solitary.py $MA/solitary/nf_sol_outer2 $MA/solitary/nw_sol \
    $OUT/ma2012_solitary.png 3
