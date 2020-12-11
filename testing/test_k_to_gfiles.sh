#!/bin/bash

# timesteps - max 16 (for 16 avail cores)
# nprocs sets the number of parallel mpi processes called - max 16 (for 16 avail cores)
nprocs=4
starttime=1000
basedir=$PWD
shot="$1"
efit_exe_pub="$2"
efit_exe="$3"
nsteps="$4"

################################################################################
################################################################################
#
################################################################################
################################################################################


if [ "$shot"x == "x" ] ; then
  echo "Requires a shot number"
  exit
fi

cd $shot/"kfiles"
echo ""
echo ""
echo "testing $efit_exe on shot $shot with/for $nsteps cores/timeslices"

################################################################################
 echo "Running public efit $shot in parallel"
################################################################################

cd public/parallel
sed -i 's|mode=.|mode=3|' efit.input

sed 's|efit_exe|'$efit_exe'|g' $basedir/efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
 echo "Running new efit $shot in parallel"
################################################################################

cd new/parallel

sed -i 's|mode=.|mode=3|' efit.input
sed 's|efit_exe|'$efit_exe'|g' $basedir/efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
echo "Running new $shot in serial"
################################################################################

cd public/serial

sed -i 's|mode=.|mode=3|' efit.input
$efit_exe 65 >> efit.log 2>&1

cd ../../


################################################################################
echo "Running public $shot in serial"
################################################################################

cd new/serial
sed -i 's|mode=.|mode=3|' efit.input

$efit_exe 65 >> efit.log 2>&1

cd ../../

################################################################################
echo " Waiting 10 seconds for job to finish"
sleep 10
echo ""
################################################################################
#wait 

np_wc=`ls -l new/parallel/g* | wc | awk '{print($1)}'`
pp_wc=`ls -l public/parallel/g* | wc | awk '{print($1)}'`
ns_wc=`ls -l new/serial/g* | wc | awk '{print($1)}'`
ps_wc=`ls -l public/serial/g* | wc | awk '{print($1)}'`

np_l=`sort new/parallel/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`
pp_l=`sort public/parallel/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`
ns_l=`sort new/serial/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`
ps_l=`sort public/serial/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`

echo "total slices of $nsteps expected"
echo " ================================="
echo "public parallel: $pp_wc files"
echo "new parallel: $np_wc files"
echo "public serial  : $ps_wc files"
echo "new serial: $ns_wc files"
echo ""

echo "Repeated timeslices (0 expected)"
echo " ================================="
echo "public parallel : $pp_l"
echo "new parallel: $np_l"
echo "public serial: $ps_l"
echo "new serial: $ns_l"
echo ""

echo "Problems and Exceptions"
echo " ================================="
echo "public parallel :"

grep "Problem" public/parallel/efit.log
grep "exception" public/parallel/efit.log
echo ""
echo "new parallel:"

grep "Problem" new/parallel/efit.log
grep "exception" new/paralel/efit.log

echo "public serial :"

grep "Problem" public/serial/efit.log
grep "exception" public/serial/efit.log
echo ""
echo "new serial:"

grep "Problem" new/serial/efit.log
grep "exception" new/serial/efit.log

cd new/serial
for t in `ls g*` ; do
  quite=`diff $t ../../public/serial/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done
cd ../..

cd new/parallel
for t in `ls g*` ; do
  quite=`diff $t ../../public/parallel/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done

echo ""
