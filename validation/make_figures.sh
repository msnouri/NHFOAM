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
#   WD   : wetting-and-drying runs (wetdry.py: $WD/final; NHWAVE runs of the
#          dam break $WD/db/nw_db, the run-up $WD/runup/nw_ru028[f] and Ting
#          & Kirby $WD/tk/nw_tk; Ting & Kirby with NHFOAM in $WD/tk)
#   BP4  : NTHMP benchmark 4 (single wave on a simple beach): profs and
#          Lab_runup.txt
#   TK94 : the measurements of Ting & Kirby (1994), Fig. 4, digitised
#          (x_TK, mean level, crest - mean, trough - mean; not distributed)
#   BRK  : breaking and longshore-current runs (Ting & Kirby: wetdry.py,
#          Hansen & Svendsen: the breakingSlope set-up, Visser case 4:
#          longshore.py; NHWAVE run of Visser case 4 in $BRK/nw_visser)
#   KW   : Ting & Kirby with the k-omega model (wetdry.tingkirby, rans=...)
#   WM   : wave-input tests (wavemaker.py)
#   UN   : unstructured meshes (unstructured.py: strip, strip_ref, cyl)
#   BAR2 : submerged bar with the boundary wavemaker (make_case.bar_case)
#   LS   : the landslide with its beach (make_case.enet_grilli_case)
#   CI   : the conical island (conical_island.py <dir>: nf_A_05 ... nw_C_10)
#   TS   : measured records on a flat tank (timeseries_input.py:
#          ts_<exact|eta|etaC|paddle>_<10|30>)
#   CB   : the composite beach (composite_beach.py: nf_<A|B|C>_<gauge4|
#          paddle|solitary>, nw_<A|B|C>)
#   MV   : the Monai valley (monai_valley.py <dir>: nf_14, nf_28, nw_28)
set -e
cd "$(dirname "$0")"
RUNS_IF=${RUNS_IF:-/home/claude/runs_if}
RUNS_CELL=${RUNS_CELL:-/home/claude/runs}
CMP=${CMP:-/home/claude/runs_cmp}
MA=${MA:-/home/claude/runs_ma}
LAB=${LAB:-/home/claude/ext/BENCHMARK_FUNWAVE}
NHW=${NHW:-/home/claude/ext/nhw_forks/NHWAVE_NHWAVE}
WD=${WD:-/home/claude/runs_wd}
BP4=${BP4:-/home/claude/ext/nthmp/nthmp-benchmark-problems/BP04-JosephZ-Single_wave_on_simple_beach}
TK94=${TK94:-/home/claude/ext/tk94/dfdm/examples/tingkirby94_spillingwaves/exp/tk94_fig4_spilling.dat}
BRK=${BRK:-/home/claude/runs_brk}
KW=${KW:-/home/claude/runs_kw}
WM=${WM:-/home/claude/runs_wm}
UN=${UN:-/home/claude/runs_un}
BAR2=${BAR2:-/home/claude/runs_bar2}
LS=${LS:-/home/claude/runs_ls}
CI=${CI:-/home/claude/runs_ci}
TS=${TS:-/home/claude/runs_ts}
CB=${CB:-/home/claude/runs_cb}
MV=${MV:-/home/claude/runs_mv}
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
python3 plot_ma2012_basin.py $MA/basin $OUT/ma2012_basin.png nfe
python3 plot_ma2012_solitary.py --old $MA/solitary/nf_sol $MA/solitary/nf_sol_ext_flt \
    $MA/solitary/nw_sol $OUT/ma2012_solitary.png 3
python3 plot_bar.py --tend 40 --title 'Submerged bar, case A ($H = 2$ cm, $T = 2.02$ s), 3 layers' \
    $LAB/car_luth_AC/lab_data/CaseA $OUT/ma2012_bar.png \
    $MA/bar/nf_bar_L3_ext_flt:NHFOAM nhwave:$MA/bar/nw_bar_L3:NHWAVE \
    $MA/bar/nf_bar_L3:"NHFOAM, first order in time"
python3 plot_enet_grilli.py $NHW/examples/landslide_benchmarks/BM2_enet_grilli/d61 \
    $OUT/ma2012_landslide.png \
    $MA/landslide/nf_eg_ext_flt:NHFOAM nhwave:$CMP/nw_eg:NHWAVE \
    $RUNS_IF/ieg_llv2:"NHFOAM before the fixes"
python3 plot_ma2012_berkhoff.py $NHW/examples/berkhoff/postprocessing $OUT/ma2012_berkhoff.png \
    $MA/berkhoff/nf_berk:"NHFOAM@32.25" nhwave:$MA/berkhoff/nw_berk:NHWAVE

# Wetting and drying (report section 6, README 6.6)
python3 plot_wetdry.py thacker $WD/final/thacker $OUT/wd_thacker.png
python3 plot_wetdry.py dambreak $OUT/wd_dambreak.png $WD/final/db:NHFOAM \
    $WD/final/db_plain:"NHFOAM, centred face depth" nhwave:$WD/db/nw_db:NHWAVE
python3 plot_runup.py $BP4/profs $BP4/Lab_runup.txt $OUT/wd_runup.png 0.0588 \
    $WD/final/ru:NHFOAM $WD/final/ru_plain:"NHFOAM, centred face depth" \
    nhwave:$WD/runup/nw_ru028:NHWAVE
python3 plot_runup.py $BP4/profs $BP4/Lab_runup.txt /tmp/wd_runup_all.png 0.0588 \
    $WD/final/ru:NHFOAM $WD/final/ru_upw:"NHFOAM, upwind face depth" \
    $WD/final/ru_plain:"NHFOAM, centred face depth" \
    $WD/final/ru_cd:"NHFOAM, \$C_d = 0.006\$" nhwave:$WD/runup/nw_ru028f:NHWAVE
cp /tmp/wd_runup_all_shoreline.png $OUT/wd_runup_shoreline.png
cp /tmp/wd_runup_all_shoreline.pdf $OUT/wd_runup_shoreline.pdf
python3 plot_tingkirby.py $TK94 $OUT/tk_spilling.png $WD/tk/tk_nf:NHFOAM \
    $WD/tk/tk_b01:"NHFOAM, \$\\beta = 0.1\$" nhwave:$WD/tk/nw_tk:NHWAVE

# Breaking with a breaking age and eddy viscosity; longshore current
TK5=${TK94%/*}/tk94_fig5_undertow_spilling.dat
python3 plot_tingkirby.py $TK94 $OUT/tk_breaking.png \
    $BRK/tk_hfa3:"NHFOAM, HFA" $BRK/tk_hybh3:"NHFOAM, HFA + eddy viscosity" \
    $BRK/tk_hyb3:"NHFOAM, HFA + eddy viscosity, vertical mixing" \
    nhwave:$WD/tk/nw_tk:NHWAVE
python3 plot_undertow.py $TK5 $OUT/tk_undertow.png \
    $BRK/tk_hfa3:"NHFOAM, HFA" $BRK/tk_hybh3:"NHFOAM, HFA + eddy viscosity" \
    $BRK/tk_hyb3:"NHFOAM, HFA + eddy viscosity, vertical mixing" \
    nhwave:$WD/tk/nw_tk:NHWAVE
python3 plot_hansen_svendsen.py $LAB/car_hansen_svendsen/postprocessing/ht_setup_061071.txt 1.667 \
    $OUT/hs_breaking.png \
    $RUNS_IF/ihs061_or:"NHFOAM, HFA" $BRK/hs_hyb3:"NHFOAM, HFA + eddy viscosity"
python3 plot_longshore.py $OUT/longshore.png --snapshot $BRK/ls_hybh \
    $BRK/ls_hybh:"NHFOAM, HFA + eddy viscosity" \
    $BRK/ls_hyb:"NHFOAM, ... with vertical mixing" $BRK/ls_hfa:"NHFOAM, HFA" \
    nhwave:$BRK/nw_visser:NHWAVE

# Wave input: boundary wavemaker, absorption, return current, irregular waves
python3 plot_wavemaker.py $WM $OUT
python3 bar_longwaves.py 1.01 $RUNS_IF/ibarC_L10:"zone" $BAR2/barC_zone_rc:"zone, return current" \
    $BAR2/barC_bc:"boundary wavemaker" nhwave:$CMP/nw_barC65:NHWAVE reef3d:$CMP/r3_barC65:REEF3D
for t in 45 50 55 60 65; do
    python3 plot_bar.py --caseC --tend $t $LAB/car_luth_AC/lab_data/CaseC /tmp/bar_bc_$t.png \
        $RUNS_IF/ibarC_L10:"zone" $BAR2/barC_zone_rc:"zone, return current" \
        $BAR2/barC_bc:"boundary wavemaker" | tail -1
done
python3 plot_bar.py --caseC --tend 65 $LAB/car_luth_AC/lab_data/CaseC $OUT/bar_bc_C.png \
    $RUNS_IF/ibarC_L10:"NHFOAM, zone without return current" \
    $BAR2/barC_zone_rc:"NHFOAM, zone with return current" \
    $BAR2/barC_bc:"NHFOAM, boundary wavemaker"

# Unstructured meshes
python3 plot_unstructured.py runup $OUT/unstructured_runup.png $WD/final/ru \
    $BP4/Lab_runup.txt $BP4/profs $UN/strip:"triangles, 2 cm" \
    $UN/strip_ref:"triangles, 3 cm to 1.5 cm"
python3 plot_unstructured.py cylinder $OUT/cylinder_diffraction.png $UN/cyl

# Landslide with its beach
python3 plot_enet_grilli.py $NHW/examples/landslide_benchmarks/BM2_enet_grilli/d61 \
    $OUT/landslide_beach.png $LS/eg_beach:"NHFOAM, beach" \
    $MA/landslide/nf_eg_ext_flt:"NHFOAM, wall at 0.2 m" nhwave:$CMP/nw_eg:NHWAVE

# Ting & Kirby with the stabilised k-omega model
TK7=${TK94%/*}/tk94_fig7_tke_spilling.dat
python3 plot_tingkirby.py $TK94 $OUT/tk_komega.png \
    $BRK/tk_hybh3:"NHFOAM, HFA + eddy viscosity (10 layers)" \
    $KW/tk_kw20:"NHFOAM, \$k\$--\$\\omega\$ (20 layers)" \
    $KW/tk_kwhfa10:"NHFOAM, HFA + \$k\$--\$\\omega\$ (10 layers)" nhwave:$WD/tk/nw_tk:NHWAVE
python3 plot_undertow.py $TK5 $OUT/tk_komega_undertow.png \
    $BRK/tk_hybh3:"HFA + eddy viscosity" $KW/tk_kw20:"\$k\$--\$\\omega\$, 20 layers" \
    $KW/tk_kwhfa10:"HFA + \$k\$--\$\\omega\$" nhwave:$WD/tk/nw_tk:NHWAVE
python3 plot_tke.py $TK7 $OUT/tk_komega_tke.png $KW/tk_kw20:"\$k\$--\$\\omega\$, 20 layers" \
    $KW/tk_kwl0:"\$\\lambda_2 = 0\$, 20 layers" $KW/tk_kwhfa10:"HFA + \$k\$--\$\\omega\$, 10 layers"

# Solitary waves on a conical island (NTHMP benchmark 6)
CIR=()
for c in A B C; do
    CIR+=("$c:nhfoam:$CI/nf_${c}_05:0.05:NHFOAM, 5 cm"
          "$c:nhfoam:$CI/nf_${c}_10:0.1:NHFOAM, 10 cm"
          "$c:nhwave:$CI/nw_${c}_10:0.1:NHWAVE, 10 cm")
done
python3 plot_conical_island.py $OUT/conical_island_gauges.png \
    $OUT/conical_island_runup.png "${CIR[@]}"

# Measured records (time-series input) and the composite beach (NTHMP
# benchmark 5)
python3 plot_timeseries_input.py $OUT/timeseries_input.png $TS
python3 plot_composite_beach.py $OUT/composite_beach_gauges.png $CB

# Run-up onto a complex 3D beach, the Monai valley (NTHMP benchmark 7)
python3 plot_monai_valley.py $OUT/monai_valley $MV/nf_14:"NHFOAM, 1.4 cm" \
    $MV/nf_28:"NHFOAM, 2.8 cm" $MV/nw_28:"NHWAVE, 2.8 cm"
