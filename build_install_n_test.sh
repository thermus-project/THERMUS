#!/bin/bash

export THERMUS=`pwd`

./cmake_n_install_thermus.sh

export THERMUS=`pwd`/install

cd install/share/doc/Thermus/tests/
root -b -q prediction.C > brut_result.txt

sed -n '/predicted values/,$p' brut_result.txt > result.txt

diff result.txt prediction_expected.txt

if [ $? -ne 0 ]; then
  echo "TEST FAILED"
  exit 1
fi

echo "Test passed"
