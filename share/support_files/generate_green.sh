exe="$1"
my_device="$2"
res="$3"
save_dir="$4"

if [ "$my_device"x == "x" ] || [ "$res"x == "x" ]; then
  my_device="all"
  res="all"
fi 

# Move necessary directories to savedir
devices=("DIII-D" "NSTX" "ITER")
echo $res
if  [ "$save_dir"x != "x" ]; then
  mkdir -p  $savedir
  for device in ${devices[@]}
  do
    cp -r $device $savedir/$device
  done
fi 

# Loop over devices/resolutions, and run efund
for device in ${devices[@]}
  do

  if [ "$my_device" == "$device" ] || [ "$my_device" == "all" ]; then 
    for dir in $device/green/*/     # list directories 
    do
      dir=${dir%*/}      # remove the trailing "/"
      echo $dir
      cd $dir
        if [ "$res" == "all" ] ; then
          for NW in 65 129 257 513
          do
            $exe $NW
          done
        else 
          $exe $NW
        fi
      cd ../../..
    done
  fi

done
