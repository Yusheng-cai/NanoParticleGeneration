#!/bin/bash

export CC=/usr/local/bin/gcc
export CXX=/usr/local/bin/g++

# specify the build directory
build_type=RELEASE
build_dir=$PWD/${build_type}/

# remove build_dir if it already exists
if [ -d $build_dir ] 
    then 
    rm -r $build_dir
fi

# make the build directory
mkdir -p $build_dir

# configure the build with cmake
cd $build_dir
cmake .. -DCMAKE_BUILD_TYPE=${build_type}

# make with 8 threads
make -j 8 
