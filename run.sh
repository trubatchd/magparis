#!/bin/bash
#set -x

mpirun -np 8 paris > tmpout # 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
