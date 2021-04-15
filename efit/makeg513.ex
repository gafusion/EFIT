#! /bin/sh
#
#  command file to make 512x512 EFIT-F90 
#
LOG=efit-f90.log
EXDIR=`dirname $0`
VER="$EXDIR/unix.ex $EXDIR/public.ex $EXDIR/standalone.ex $EXDIR/new.ex"
# cd ../u   if in x directory
if [ -z "$1" ]
then
	echo "Usage: $0 <source file>"
	exit 1
fi

case $1 in

modules-efitx.f90)
      ex  modules-efitx.f90   << EOF > $LOG
g/\!higg/s//     /g
g/\!513g/s//     /g
wq! modules-efit.f90
EOF
;;
#
*)
	echo "$1 unknown"
	exit 1
	;;
esac
