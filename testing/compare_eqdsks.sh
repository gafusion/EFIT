#!/bin/bash
shot="$1"

if [ "$shot"x == "x" ] ; then
  echo "Requires a shot number"
  exit
fi
module purge
module load omfit/unstable
python3 compare_eqdsks.py $shot
