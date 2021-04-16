#!/bin/bash

# choose your executable (full path)
#efit_exe="efitd90"
#efit_exe="/home/kostukm/gitrepos/efit-bugtesting/efitbuild/efitdmpi90"
efit_exe="/home/kostukm/gitrepos/efit/efitbuild/efitd90"
# nsteps sets the timesteps 
nsteps=1
# nprocs sets the number of parallel mpi processes called - max 16 (for 16 avail cores)
nprocs=1
starttime=1380

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

if [ -d "$shot" ] ; then 
  cd $shot
else
  mkdir $shot
  cd $shot
fi
###
# #
###
echo ""
echo "executable : `which $efit_exe`"
echo "shot       : $shot"
echo "timeslices : $nsteps"
echo "mpi_procs  : $nprocs"
echo "start_time : $starttime"
echo ""
module list
echo ""

rm -rf serial
rm -rf parallel

mkdir serial
mkdir parallel

################################################################################
 echo "Running Parallel $shot"
################################################################################

cd parallel
date >> efit.log 

sed 's|nsteps|'$nsteps'|g' ../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../efit_snap* ./

sed 's|efit_exe|'$efit_exe'|g' ../../efit.sbatch >> efit.sbatch
sed -i 's|nprocs|'$nprocs'|g' efit.sbatch
sbatch ./efit.sbatch

cd ../

################################################################################
echo "Running Serial $shot"
################################################################################

cd serial
date >> efit.log

module list >> efit.log 2>&1
sed 's|nsteps|'$nsteps'|g' ../../efit.input >> efit.input
sed -i 's|nshot|'$shot'|g' efit.input
sed -i 's|start_time|'$starttime'|g' efit.input
cp ../../efit_snap* ./

$efit_exe 65 >> efit.log 2>&1

cd ../

################################################################################
echo " Waiting 10 seconds for job to finish"
sleep 10
################################################################################

s_wc=`ls -l serial/g* | wc | awk '{print($1)}'`
p_wc=`ls -l parallel/g* | wc | awk '{print($1)}'`

s_l=`sort serial/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`
p_l=`sort parallel/efit.log | grep "time" | uniq -d | wc | awk '{print($1)}'`

#sort parallel/efit.log 
#sort parallel/efit.log | grep time | uniq -d

#echo "Serial: $s_wc files   Parallel: $p_wc files"
echo ""
echo "total slices of $nsteps expected"
echo " ================================="
echo "Serial  : $s_wc files"
echo "Parallel: $p_wc files"
echo ""

echo "Repeated timeslices (0 expected)"
echo " ================================="
echo "Serial  : $s_l"
echo "Parallel: $p_l (>0 indicates likely serial running)"
echo ""

echo "Problems and Exceptions"
echo " ================================="
echo "Serial  :"
#ls -l serial/g*
grep "Problem" serial/efit.log
grep "exception" serial/efit.log
#grep "Error" serial/efit.log ## every fit error, not necessarily an EFIT runtime error
echo ""
echo "Parallel:"
#ls -l parallel/g*
grep "Problem" parallel/efit.log
grep "exception" parallel/efit.log
#grep "Error" parallel/efit.log ## every fit error, not necessarily an EFIT runtime error

echo ""
echo "Diff Comparisons (0 == identical, 1 == differences, 2 == file DNE)"
echo " ================================="
cd serial
for t in `ls g*` ; do
  quite=`diff $t ../parallel/$t 2>&1`
  error=$?
  echo "$t diff returns $error"
done
echo ""

#s_pb=
#diff serial/g"$shot".00100
