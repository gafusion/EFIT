#!/bin/bash

# choose your executable (full path)

efit_exe_pub="efitd90"
efit_exe="/home/mcclenaghanj/efit_testing/efitbuild/efitdmpipgf90"

# timesteps - max 16 (for 16 avail cores)
# nprocs sets the number of parallel mpi processes called - max 16 (for 16 avail cores)
nprocs=2
starttime=3000
no_snap_ext=true



################################################################################
################################################################################
#
################################################################################
################################################################################
shot="$1"

if [ "$shot"x == "x" ] ; then
  echo "Requires a shot number"
  exit
fi

if [ -d "kfiles_""$shot" ] ; then
  rm -r "kfiles_"$shot 
fi
mkdir "kfiles_"$shot
cd "kfiles_"$shot

###
# #
###
echo ""
echo ""
echo "testing $efit_exe on shot $shot with/for $nsteps cores/timeslices"
echo ""
module list
echo ""

mkdir -p new/parallel
mkdir -p public/parallel
mkdir -p new/serial
mkdir -p public/serial

################################################################################
 echo "Running public efit $shot in parallel"
################################################################################

cd public/parallel
date >> efit.log 

sed 's|nsteps|'$nsteps'|g' ../../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../../efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
## testing no snap ext
if $no_snap_ext ; then
   sed -i 's|jta_f||' efit.input
   mv efit_snap* ./efit_snap.dat
fi

sed 's|efit_exe|'$efit_exe'|g' ../../../efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
 echo "Running new efit $shot in parallel"
################################################################################

cd new/parallel
date >> efit.log 

sed 's|nsteps|'$nsteps'|g' ../../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../../efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
## testing no snap ext
if $no_snap_ext ; then
   sed -i 's|jta_f||' efit.input
   mv efit_snap* ./efit_snap.dat
fi

sed 's|efit_exe|'$efit_exe'|g' ../../../efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../../


################################################################################
echo "Running new $shot in serial"
################################################################################

cd public/serial
date >> efit.log

module list >> efit.log 2>&1
sed 's|nsteps|'$nsteps'|g' ../../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../../efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
## testing no snap ext
if $no_snap_ext ; then
   sed -i 's|jta_f||' efit.input
   mv efit_snap* ./efit_snap.dat
fi

$efit_exe 65 >> efit.log 2>&1

cd ../../

################################################################################
echo "Running public $shot in serial"
################################################################################

cd new/serial
date >> efit.log

module list >> efit.log 2>&1
sed 's|nsteps|'$nsteps'|g' ../../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../../efit_snap* ./

sed -i 's|mode=.|mode=5|' efit.input
## testing no snap ext
if $no_snap_ext ; then
   sed -i 's|jta_f||' efit.input
   mv efit_snap* ./efit_snap.dat
fi

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

echo ""
echo "Diff Comparisons (0 == identical, 1 == differences, 2 == file DNE)"
echo " ================================="
cd serial
for t in `ls k*` ; do
  quite=`diff $t ../parallel/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done
echo ""
