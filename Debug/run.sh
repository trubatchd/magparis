#!/bin/bash
#set -x

mpirun -np 1 magparis > tmpout # 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
