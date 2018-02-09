#!/bin/bash -e
# checked on 30/01/2018 with root-v6-12-04 (B.H.)
echo $THERMUS
export BASEDIR=$(pwd)
#
if [ ! $THERMUS ]
  then echo 'first make sure $THERMUS is set to the THERMUS code directory'
else
#
cd $THERMUS
rm -rf build
mkdir build
cd build
cmake -Wdev -debug-output -DCMAKE_VERBOSE_MAKEFILE=ON ..
cmake --build . -- -j10
cd $BASEDIR
#
fi
