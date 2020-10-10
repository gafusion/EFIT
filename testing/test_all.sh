#!/bin/bash

for shot in 160606 172461 182943
do
rm -rf "magnetics_"$shot
rm -rf "kfiles_"$shot
module purge
module load efit
./test_magnetics.sh $shot 
./compare_eqdsks.sh "magnetics_"$shot

./test_kfiles.sh $shot

done
