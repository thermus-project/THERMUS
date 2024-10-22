#!/bin/bash -x
# This file is not distributed, it is used to test THERMUS
# before committing changes
TESTDIR=`dirname -- "$( readlink -f -- "$0"; )"`
BASEDIR=`dirname $TESTDIR`

EXTRA_OPTIONS="$@"

SCRIPTDIR=$BASEDIR/scripts

# Remove old builds and create build directory
rm -rf $BASEDIR/build $BASEDIR/local
rm -f $BASEDIR/brut_result.txt $BASEDIR/result.txt

$SCRIPTDIR/inplace_build.sh $EXTRA_OPTIONS

$TESTDIR/inplace_test.sh $BASEDIR/run_thermus

