#!/bin/sh

cd RGBDPlaneDetection/include/MRF2.2
make
cd ../../..
mkdir build
cd build
cmake ../RGBDPlaneDetection
make
cd ..