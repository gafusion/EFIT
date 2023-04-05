exe="$1"
my_device="$2"
res="$3"
save_dir="$4"

if [ "$my_device"x == "x" ] || [ "$res"x == "x" ]; then
  my_device="all"
  res="all"
fi 

echo $res
if  [ "$save_dir"x != "x" ]; then
  mkdir -p  $savedir
  for device in "DIII-D" "NSTX" "ITER"
  do
    cp -r $device $savedir/$devices
  done
fi 


for device in "DIII-D" "NSTX" "ITER"
  do

  if [ "$my_device" == "$device" ] || [ "$my_device" == "all" ]; then 
    for dir in $device/green/*/     # list directories in the form "/tmp/dirname/"
    do
      dir=${dir%*/}      # remove the trailing "/"
      echo $dir
      cd $dir
        if [ "$res" == "all" ] ; then
          for NW in 65 129 257 513
          do
            sed -i 's|NW = .*|NW = '$NW'|g' mhdin.dat
            $exe
          done
        else 
          sed -i 's|NW = .*|NW = '$NW'|g' mhdin.dat
          $exe
        fi
      cd ../../..
    done
  fi

done
