#!/bin/bash
# Period error with surface-refined sigma layers (blockMesh simpleGrading in z:
# ratio of the top to the bottom layer thickness). Needs OpenFOAM + NHFoam.
here=$(dirname $(readlink -f $0))
out=./grading.csv
echo "kh,layers,grading,T_lin,T_nh,ratio,amp_ratio" > $out
run_one() {
  kh=$1; nl=$2; gr=$3; name=grade_kh${kh}_L${nl}_g${gr}
  rm -rf $name
  python3 - << PY
import sys; sys.path.insert(0,'$here')
import make_case
r = make_case.standing_case('$name', $kh, $nl, Nx=64, stepsPerPeriod=200, nPeriods=5, zgrade=$gr)
open('$name/T_lin','w').write(repr(r['T']))
PY
  (cd $name && blockMesh > log.blockMesh 2>&1 && nhSetWaves > log.nhSetWaves 2>&1 && NHFoam > log.NHFoam 2>&1)
  res=$(python3 $here/analyse_standing.py $name $(cat $name/T_lin))
  Tn=$(echo "$res" | sed -n 's/.*period NHFoam = \([0-9.]*\).*/\1/p'); Tl=$(cat $name/T_lin)
  rat=$(echo "$res" | sed -n 's/.*ratio = \([0-9.]*\)$/\1/p' | head -1)
  amp=$(echo "$res" | sed -n 's/.*periods: \([0-9.]*\)/\1/p')
  echo "$kh,$nl,$gr,$Tl,$Tn,$rat,$amp" >> $out
  rm -rf $name
}
for kh in 1 2 3 5; do
  for nl in 3 5; do
    for gr in 1 0.5 0.25; do
      run_one $kh $nl $gr &
      while [ $(jobs -r | wc -l) -ge 2 ]; do sleep 0.2; done
    done
  done
done
wait
sort -t, -k1,1g -k2,2n -k3,3gr $out
