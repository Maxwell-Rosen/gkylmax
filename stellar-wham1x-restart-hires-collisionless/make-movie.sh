#!/bin/bash
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: `basename $0`"
  echo "-r frame rate, r frames per second. defaults to 2"
  echo "-i input files, e.g., foo_%04d.png, %04d ==> 4 digit file extension. No default value"
  echo "-o output name, e.g., foo.mp4. No default value"
  exit 0
fi

#Default rate
rate=2

while [ -n "$1" ]
 
# while loop starts
 
do
 
case "$1" in
 
-r) rate="$2";;
 
-i) input="$2";;
  
-o) output="$2";;
 
--) shift
 
break ;;
  
esac
 
shift
 
done

ffmpeg -f image2 -r $rate -i $input -c:v libx264 -pix_fmt yuv420p -y  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $output