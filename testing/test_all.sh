#!/bin/bash


efit_exe_pub="efitd90"
efit_exe="/home/mcclenaghanj/efit_testing/efitbuild/efitdmpipgf90"
rm -r logs
mkdir logs

for shot in 163303 172461 182943

do
module purge
module load efit
rm -r $shot

echo "Testing EFIT01 for shot  $shot"
./test_snap.sh $shot "EFIT01" $efit_exe_pub $efit_exe > logs/test_EFIT01_$shot.log
./compare_eqdsks.sh "$shot/EFIT01" >> logs/test_EFIT01_$shot.log
grep "Maximum difference in" logs/test_EFIT01_$shot.log
echo ""

echo "Testing EFIT02 for shot $shot"
./test_snap.sh $shot "EFIT02" $efit_exe_pub $efit_exe > logs/test_EFIT02_$shot.log
./compare_eqdsks.sh "$shot/EFIT02" >> logs/test_EFIT02_$shot.log
grep "Maximum difference in" logs/test_EFIT02_$shot.log
echo ""

echo " Testing kfiles generation and gfile produced from kfiles for shot  $shot"
./test_kfiles.sh $shot 'EFIT01' $efit_exe_pub $efit_exe > logs/test_kfiles_$shot.log
./test_k_to_gfiles.sh $shot $efit_exe_pub $efit_exe > logs/test_k_to_gfiles_$shot.log
./compare_eqdsks.sh "$shot/kfiles" >> logs/test_k_to_gfiles_$shot.log
grep "Maximum difference in" logs/test_k_to_gfiles_$shot.log
echo ""

done

echo " Testing gfile generated from kineticEFIT kfile"
./test_kineticEFIT.sh  $efit_exe_pub $efit_exe > logs/test_kineticEFIT.log

