#!/bin/bash -e
# checked on 20/07/2021 with root-v6-24-02 (B.H.)
nproc=`grep -c "^processor" /proc/cpuinfo`

echo $THERMUS
export BASEDIR=$(pwd)
#
if [[ (`export | egrep "THERMUS\=" | wc -l` == 0 ) || ( ! $THERMUS ) ]]
  then echo 'first make sure the shell variable ${THERMUS} is set to the THERMUS code directory and properly exported'
  exit -1
fi

#
cd $THERMUS
rm -rf build install
mkdir build
cd build
#cmake -Wdev -debug-output -DCMAKE_VERBOSE_MAKEFILE=ON ..
cmake -Wdev --debug-output -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX=$BASEDIR/install ..
cmake --build . --parallel $nproc
cmake --install .

# FIXME: temporary for check
ls -lR $BASEDIR/install
