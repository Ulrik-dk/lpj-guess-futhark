#!/bin/bash

# Need to be adapted to your HPC environment:
# ml intel/2018b
# ml netCDF/4.6.1-intel-2018b
# ml netCDF-C++4/4.3.0-intel-2018b
# ml CMake/3.12.1-GCCcore-7.3.0

#CC=icc
#CXX=icc

CC=gcc
CXX=gcc

cd source/build
rm CMakeCache.txt

# The line below is important, otherwise you probably won't be able to read CF files.
# Find the location of the netCDF library on your HPC, probably based on the same compiler
# as the modules loaded at the top of this script (eg. intel-2018b in our case).
# When you run this compile script, in the beginning it should mention "Found NetCDF"
# CMAKE_PREFIX_PATH=/apps/gent/CO7/skylake-ib/software/netCDF/4.6.1-intel-2018b/lib64/
CMAKE_PREFIX_PATH=/usr/lib/x86_64-linux-gnu/

cmake ..
make

# Comment this out if you don't want to overwrite the binary
# in the example folder automatically...
cp guess ../../example/

