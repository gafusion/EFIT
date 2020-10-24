#!/bin/bash

# timesteps - max 16 (for 16 avail cores)
# nprocs sets the number of parallel mpi processes called - max 16 (for 16 avail cores)
nsteps=1
nprocs=1
starttime=1400
shot=164409
basedir=$PWD
efit_exe_pub="$1"
efit_exe="$2"

################################################################################
################################################################################
#
################################################################################
################################################################################


if [ -d "$shot/""kineticEFIT" ] ; then
  rm -r $shot/"kineticEFIT" 
fi
mkdir -p $shot/"kineticEFIT"
cd $shot/"kineticEFIT"

mkdir -p new/serial
mkdir -p public/serial

echo ""
echo ""
echo "testing kinetic EFIT $efit_exe on shot $shot with/for $nsteps cores/timeslices"


################################################################################
echo "Running new $shot in serial"
################################################################################

cd public/serial

cp $basedir/k$shot.0$starttime .
sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
sed -i 's|mode=.|mode=3|' efit.input
$efit_exe 129 >> efit.log 2>&1

cd ../../


################################################################################
echo "Running public $shot in serial"
################################################################################

cd new/serial
cp $basedir/k$shot.0$starttime .
sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
sed -i 's|mode=.|mode=3|' efit.input

$efit_exe 129 >> efit.log 2>&1

cd ../../

ns_wc=`ls -l new/serial/g* | wc | awk '{print($1)}'`
ps_wc=`ls -l public/serial/g* | wc | awk '{print($1)}'`

ns_l=`sort new/serial/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`
ps_l=`sort public/serial/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`

echo "total slices of $nsteps expected"
echo " ================================="
echo "public serial  : $ps_wc files"
echo "new serial: $ns_wc files"
echo ""

echo "Repeated timeslices (0 expected)"
echo " ================================="
echo "public serial: $ps_l"
echo "new serial: $ns_l"
echo ""

echo "Problems and Exceptions"
echo " ================================="

echo "public serial :"

grep "Problem" public/serial/efit.log
grep "exception" public/serial/efit.log
echo ""
echo "new serial:"

grep "Problem" new/serial/efit.log
grep "exception" new/serial/efit.log
echo ""
cd new/serial
for t in `ls g*` ; do
  quite=`diff $t ../../public/serial/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done
echo ""
