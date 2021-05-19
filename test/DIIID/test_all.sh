##!/bin/bash


#module purge
#module load env/gcc.9.2
export efit="/home/mcclenaghanj/efit_ai/efit/efit/efit"
export link_efit="/fusion/projects/theory/mcclenaghanj/efit_support_files/d3d/"

for dir in efit01 efit02 rfile kineticEFIT

do 
echo "Testing" $dir
cd $dir
./run_command.sh 
cd ..
echo "Done testing" $dir

done
