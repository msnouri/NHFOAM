#!/bin/bash
# requires the OpenFOAM environment and the NHFOAM executables
# PRESSURE=cell|interface selects the pressure layout (default interface)
out=./dispersion.csv
echo "kh,layers,T_lin,T_nh,ratio,amp_ratio" > $out
run_one() {
  kh=$1; nl=$2; name=disp_kh${kh}_L${nl}
  rm -rf $name
  python3 - << PY
import sys; sys.path.insert(0,'$(dirname $(readlink -f $0))')
import make_case, os
make_case.PRESSURE = os.environ.get('PRESSURE', 'interface')
r = make_case.standing_case('$name', $kh, $nl, Nx=64, stepsPerPeriod=200, nPeriods=5)
open('$name/T_lin','w').write(repr(r['T']))
PY
  (cd $name && blockMesh > log.blockMesh 2>&1 && nhSetWaves > log.nhSetWaves 2>&1 && NHFOAM > log.NHFOAM 2>&1)
  res=$(python3 $(dirname $(readlink -f $0))/analyse_standing.py $name $(cat $name/T_lin))
  Tn=$(echo "$res" | sed -n 's/.*period NHFOAM = \([0-9.]*\).*/\1/p'); Tl=$(cat $name/T_lin)
  rat=$(echo "$res" | sed -n 's/.*ratio = \([0-9.]*\)$/\1/p' | head -1)
  amp=$(echo "$res" | sed -n 's/.*periods: \([0-9.]*\)/\1/p')
  echo "$kh,$nl,$Tl,$Tn,$rat,$amp" >> $out
  rm -rf $name/[1-9]* 
}
for kh in 0.5 1 2 3 5 7; do
  for nl in 1 2 3 5 10 20; do
    run_one $kh $nl &
    while [ $(jobs -r | wc -l) -ge 2 ]; do sleep 0.2; done
  done
done
wait
sort -t, -k1,1g -k2,2n $out | column -t -s,
