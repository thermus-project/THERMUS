#!/bin/bash -x
# Inplace test for THERMUS
# As we don't know where run_thermus is located, we need to pass it as an argument
RUN_THERMUS=$1

# Location of test files
TESTDIR=`dirname -- "$( readlink -f -- "$0"; )"`

# Remove old tests resulst
rm -f $TESTDIR/results/*.txt

# Run test
# $THERMUS is not yet defined
# So the argument is single-quoted to avoid expansion
$RUN_THERMUS '$THERMUS/share/doc/Thermus/tests/all_predictions.C -b -q'

# If there are no results files, test failed
if compgen -G "$TESTDIR/results/*.txt" > /dev/null; then
    echo "Some results files were found"
else
    echo "TEST FAILED: no results files were found"
    exit 1
fi

# For each file in the result directory, test if it's same as in the expected directory
for ff in `ls $TESTDIR/results/*.txt`; do
  f=`basename $ff`
  diff $TESTDIR/results/$f $TESTDIR/expected_results/$f
  if [ $? -ne 0 ]; then
    echo "TEST FAILED for file $f"
    exit 1
  fi
done

echo "Test passed"
