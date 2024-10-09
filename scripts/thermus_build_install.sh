#!/bin/bash -e
# Build and install THERMUS
# If no argument given, install will be done in $BASEDIR/local
SCRIPTDIR=`dirname -- "$( readlink -f -- "$0"; )"`
BASEDIR=`dirname $SCRIPTDIR`

INSTALL_DIR=${1:-$BASEDIR/local}

EXTRA_OPTIONS="${@:2}"

# Get number of processors for parallel build
nproc=`${SCRIPTDIR}/ncore.sh -l`

# Remove old builds and create build directory
cd $BASEDIR
rm -rf build local
mkdir build

# Invoque cmake in build directory
cmake -Wdev --debug-output -DCMAKE_VERBOSE_MAKEFILE=ON -B build -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR $EXTRA_OPTIONS
cmake --build build --parallel $nproc
cmake --install build

# List installed files
ls -lR $INSTALL_DIR
