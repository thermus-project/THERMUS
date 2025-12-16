#!/bin/bash -x
# This file is not distributed, it is used to test THERMUS
# before committing changes
MASTERSCRIPTDIR=`dirname -- "$( readlink -f -- "$0"; )"`
MASTERDIR=`dirname $MASTERSCRIPTDIR`
WORKDIR=$MASTERDIR/install_test

rm -rf $WORKDIR

THERMUSGIT=$MASTERDIR
THERMUSBRANCH="dev4ywork6"

mkdir -p $WORKDIR/src
mkdir -p $WORKDIR/local
mkdir -p $WORKDIR/tests
cd $WORKDIR/src
git clone $THERMUSGIT

cd THERMUS
git checkout $THERMUSBRANCH

./scripts/thermus_build_install.sh $WORKDIR/local

cd $WORKDIR/tests
rm -rf $WORKDIR/src

export PATH=$WORKDIR/local/bin:$PATH

$WORKDIR/local/share/doc/Thermus/tests/test_on_installed.sh

