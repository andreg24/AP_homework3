#!/bin/bash

# builds dynamic libraries???

set -x # print commands executed

BOOST_ROOT=~/Desktop/uni/advanced_programming/homeworks/extlib # $1
EIGEN_PATH=/usr/local/include/eigen3 # $2
# PYBIND_PATH=$3 # example: ... already specified in CmakeLists.txt

# execute cmake command

if [ -d "./build" ]; then # checks if folder exist
	rm -r build
fi
mkdir build
cd build
cmake -DBOOST_ROOT=${BOOST_ROOT} -DEIGEN_PATH=${EIGEN_PATH} ..
make
cd ..
