#!/bin/bash

# choose your executable (full path)



# timesteps - max 16 (for 16 avail cores)
# nprocs sets the number of parallel mpi processes called - max 16 (for 16 avail cores)
nsteps=4
nprocs=4
starttime=1000
basedir=$PWD
shot="$1"
snapext="$2"
efit_exe_pub="$3"
efit_exe="$4"

echo $basedir
################################################################################
################################################################################
#
################################################################################
################################################################################


if [ "$shot"x == "x" ] ; then
  echo "Requires a shot number"
  exit
fi

if [ -d "$shot/""kfiles" ] ; then
  rm -r $shot/"kfiles"
fi

mkdir -p $shot/"kfiles"
cd $shot/"kfiles"

###
# #
###
echo ""
echo ""
echo "testing $efit_exe on shot $shot with/for $nsteps cores/timeslices"

mkdir -p new/parallel
mkdir -p public/parallel
mkdir -p new/serial
mkdir -p public/serial

################################################################################
 echo "Running public efit $shot in parallel"
################################################################################

cd public/parallel
date >> efit.log 

sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp $basedir/efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
sed -i 's|jta_f|'$snapext'|g' efit.input

sed 's|efit_exe|'$efit_exe'|g' $basedir/efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
 echo "Running new efit $shot in parallel"
################################################################################

cd new/parallel
date >> efit.log 

sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp $basedir/efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
sed -i 's|jta_f|'$snapext'|g' efit.input

sed 's|efit_exe|'$efit_exe'|g' $basedir/efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
echo "Running new $shot in serial"
################################################################################

cd public/serial
date >> efit.log

module list >> efit.log 2>&1
sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp $basedir/efit_snap* ./
sed -i 's|mode=.|mode=5|' efit.input
sed -i 's|jta_f|'$snapext'|g' efit.input

$efit_exe 65 >> efit.log 2>&1

cd ../../

################################################################################
echo "Running public $shot in serial"
################################################################################

cd new/serial
date >> efit.log

module list >> efit.log 2>&1
sed 's|nsteps|'$nsteps'|g' $basedir/efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp $basedir/efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
sed -i 's|jta_f|'$snapext'|g' efit.input

$efit_exe 65 >> efit.log 2>&1

cd ../../

################################################################################
echo " Waiting 10 seconds for job to finish"
sleep 10
echo ""
################################################################################
#wait 

np_wc=`ls -l new/parallel/k* | wc | awk '{print($1)}'`
pp_wc=`ls -l public/parallel/k* | wc | awk '{print($1)}'`
ns_wc=`ls -l new/serial/k* | wc | awk '{print($1)}'`
ps_wc=`ls -l public/serial/k* | wc | awk '{print($1)}'`

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
grep "exception" new/parallel/efit.log

echo "public serial :"

grep "Problem" public/serial/efit.log
grep "exception" public/serial/efit.log
echo ""
echo "new serial:"

grep "Problem" new/serial/efit.log
grep "exception" new/serial/efit.log

echo ""
echo "Diff Comparisons (0 == identical, 1 == differences, 2 == file DNE)"
echo " ================================="

cd new/serial
for t in `ls k*` ; do
  quite=`diff $t ../../public/serial/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done
cd ../..

cd new/parallel
for t in `ls k*` ; do
  quite=`diff $t ../../public/parallel/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done

echo ""
