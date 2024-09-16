#!/bin/bash

export THERMUS=`pwd`
./cmakethermus.sh
root -b -q test/prediction.C > brut_result.txt

sed -n '/predicted values/,$p' brut_result.txt > result.txt

diff result.txt test/prediction_expected.txt

if [ $? -ne 0 ]; then
  echo "TEST FAILED"
  exit 1
fi

echo "Test passed"
