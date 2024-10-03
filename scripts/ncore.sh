#!/usr/bin/env sh

# Taken from: https://gist.github.com/Leedehai/b7ebf018b28eaadceed957072690245e
# Copyright: see README and LICENSE under the project root directory.
# Author: @Leedehai
#
# File: ncore.sh
# ---------------------------
# This is a handy script to get CPU core number on Linux and Darwin (macOS).

# Logical CPU cores
LCORE_COUNT=$([ $(uname) = 'Darwin' ] && sysctl -n hw.logicalcpu_max || lscpu -p | egrep -v '^#' | wc -l)
# Physical CPU cores
PCORE_COUNT=$([ $(uname) = 'Darwin' ] && sysctl -n hw.physicalcpu_max || lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)

if [ $# -ne 1 ]; then
	echo "One option needed. --help for more information."
	exit 0
fi

OPTION=$1

if [ $OPTION = "-h" ] || [ $OPTION = "--help" ]; then
	echo "Get CPU core number on Linux and Darwin (macOS)"
	echo "  -h | --help        show this help message"
	echo "  -l | --logical     show logical CPU count"
	echo "  -p | --physical    show physical CPU count"
	echo "  -a | --all         show both above"
elif [ $OPTION = "-l" ] || [ $OPTION = "--logical" ]; then
	echo $LCORE_COUNT
elif [ $OPTION = "-p" ] || [ $OPTION = "--physical" ]; then
	echo $PCORE_COUNT
elif [ $OPTION = "-a" ] || [ $OPTION = "--all" ]; then
	echo "$LCORE_COUNT $PCORE_COUNT"
else
	echo "Wrong option. --help for more information."
	exit 1
fi

