#!/bin/sh
# Function to show colored text on screen
cecho() {
  local code="\033["
  case "$1" in
    black  | bk) color="${code}0;30m";;
    red    |  r) color="${code}1;31m";;
    green  |  g) color="${code}1;32m";;
    yellow |  y) color="${code}1;33m";;
    blue   |  b) color="${code}1;34m";;
    purple |  p) color="${code}1;35m";;
    cyan   |  c) color="${code}1;36m";;
    gray   | gr) color="${code}0;37m";;
    *) local text="$1"
  esac
  [ -z "$text" ] && local text="$color$2${code}0m"
  echo -e "$text"
}

# RGBD="/d/3drecon/data/3dlite/apt/apt-rgbd-render/"
# RGBD="/d/3drecon/data/bundlefusion/copyroom/copyroom/"
# OUTPUT="/d/3drecon/data/bundlefusion/copyroom/copyroom-plane-temp/"
ORIGINALRGBD="/d/3drecon/data/bundlefusion/copyroom/copyroom/"
RGBD="/d/3drecon/data/bundlefusion/office0/office0_temp/office0_temp-original-clusters800-inner0.001s-border0.1-clean/"
OUTPUT=$RGBD
if [ ! -d "$OUTPUT" ]; then
  mkdir $OUTPUT
fi
code="/c/dev/RGBDPlaneDetection/x64/Release/RGBDPlaneDetection.exe"
cecho g "input rgbd folder: $RGBD"
cecho g "output plane folder: $OUTPUT"
cecho g "RGBDPlaneDetection code: $code"
# for file in `ls $RGBD*-depth.png | sort`; do
#     echo $file >> list.txt
#     python $genply 
# done
num=200
for i in  $(seq -f "%06g" 0 $num)
do
	cecho g "Reading frame $i"
  # $code ${RGBD}frame-$i-color.png ${RGBD}frame-$i-depth.png $OUTPUT
  # $code -o ${RGBD}frame-$i.color.jpg ${RGBD}frame-$i.depth.png $OUTPUT
  $code ${RGBD}frame-$i.rendered-color.jpg ${RGBD}frame-$i.rendered-depth.png $OUTPUT
	# $code -o ${RGBD}frame-$i-color.png ${RGBD}frame-$i-depth.png
  cecho g "Copying pose file of frame $i"
  # cp ${ORIGINALRGBD}frame-$i.pose.txt ${OUTPUT}frame-$i-pose.txt
done


