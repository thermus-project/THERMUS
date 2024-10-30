#!/bin/bash -x
# Inplace test for THERMUS
# As we don't know where run_thermus is located, we need to pass it as an argument
RUN_THERMUS=$1

# Location of test files
TESTDIR=`dirname -- "$( readlink -f -- "$0"; )"`

# $THERMUS is not yet defined
# So the argument is single-quoted to avoid expansion
$RUN_THERMUS '$THERMUS/share/doc/Thermus/tests/prediction.C -b -q' > brut_result.txt

sed -n '/predicted values/,$p' brut_result.txt > result.txt

diff result.txt $TESTDIR/prediction_expected.txt

if [ $? -ne 0 ]; then
  echo "TEST FAILED"
  exit 1
fi

echo "Test passed"
