#!/bin/bash

export THERMUS=`pwd`

# Doing this here, as then the directory will exists even
# when build fails
rm -rf install
mkdir install

./cmake_n_install_thermus.sh
rm -rf brut_result.txt result.txt

cd install
root -b -q test/prediction.C > brut_result.txt

sed -n '/predicted values/,$p' brut_result.txt > result.txt

diff result.txt test/prediction_expected.txt

if [ $? -ne 0 ]; then
  echo "TEST FAILED"
  exit 1
fi

echo "Test passed"
