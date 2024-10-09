#!/bin/bash -e
SCRIPTDIR=`dirname -- "$( readlink -f -- "$0"; )"`
BASEDIR=`dirname $SCRIPTDIR`

EXTRA_OPTIONS="$@"

$SCRIPTDIR/thermus_build_install.sh $BASEDIR/local $EXTRA_OPTIONS

# Create symlink to run_thermus for easier use
ln -f -s $BASEDIR/local/bin/run_thermus $BASEDIR/run_thermus
