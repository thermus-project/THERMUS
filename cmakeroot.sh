#!/bin/bash -e
# # checked 07/03/2018 with root-v6-12-06 (B.H.)
export BASEDIR=$(pwd)
export ROOTDIR=$BASEDIR/root
cd $ROOTDIR
export ROOTVERSION=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
cd $BASEDIR
export ROOTBUILDDIR="$BASEDIR/root-build-$ROOTVERSION"
export ROOTINSTDIR="$BASEDIR/root-inst-$ROOTVERSION"
rm -rf $ROOTBUILDDIR $ROOTINSTDIR
mkdir $ROOTBUILDDIR $ROOTINSTDIR
cd $ROOTBUILDDIR
cmake $ROOTDIR -DCMAKE_INSTALL_PREFIX=$ROOTINSTDIR -Dgsl_shared:BOOL=ON
cmake --build . --target install -- -j10
cd $BASEDIR
source $ROOTINSTDIR/bin/thisroot.sh
