#!/bin/bash -x
# Installed test for THERMUS
# As we don't know where run_thermus is located, we need to pass it as an argument
`run_thermus --getenv`

# Location of test files
TESTDIR=$THERMUS/share/doc/Thermus/tests

RESULTDIR=`pwd`/tests/results
mkdir -p $RESULTDIR
# Remove old tests resulst
rm -f $RESULTDIR/*.txt

# Create test results directory
mkdir -p $RESULTDIR

# So the argument is single-quoted to avoid expansion
run_thermus $THERMUS/share/doc/Thermus/tests/all_predictions.C -b -q > $RESULTDIR/brut_result.txt
# Brut results have to be massaged to extract the final part we compare:
sed -n '/predicted values/,$p' $RESULTDIR/brut_result.txt > $RESULTDIR/result.txt
rm $RESULTDIR/brut_result.txt # Needed as we count output files

# LHC5020 test
run_thermus '$THERMUS/share/doc/Thermus/tests/lhc5020_fit_charm.C -b -q'
# If there are no results files, test failed
if compgen -G "$RESULTDIR/*.txt" > /dev/null; then
    echo "Some results files were found"
else
    echo "TEST FAILED: no results files were found"
    exit 1
fi

# Number of test files
nbf=`compgen -G "$RESULTDIR/*.txt" | wc -l`

# Number of expected files
nbe=`compgen -G "$THERMUS/share/doc/Thermus/tests/expected_results/*.txt" | wc -l`

# All comparisons will be performed, even if some failed, this is to ensure better debugging
NBF_FAILED=0
if [ $nbf -ne $nbe ]; then
  echo "TEST FAILED: number of expected files does not match number of results files"
  NBF_FAILED=1
fi

RES_FAILED=0
# For each file in the result directory, test if it's same as in the expected directory
for ff in `ls $RESULTDIR/*.txt`; do
  f=`basename $ff`
  diff $RESULTDIR/$f $THERMUS/share/doc/Thermus/tests/expected_results/$f
  if [ $? -ne 0 ]; then
    echo "TEST FAILED for file $f"
    RES_FAILED=1
  fi
done

if [ $NBF_FAILED -eq 1 ] || [ $RES_FAILED -eq 1 ]; then
  # Print some debugging information and exit with an error
  env
  root --version
  cc --version
  cpp --version

  echo "TEST FAILED"
  exit 1
fi

echo "Test passed"
