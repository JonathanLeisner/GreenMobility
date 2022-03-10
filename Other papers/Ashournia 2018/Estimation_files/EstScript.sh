#!/bin/sh
#--------------------------------------------------------
# By Damoun Ashournia, University of Oxford
#--------------------------------------------------------
# Bash script to compile and run the estimation module.
#--------------------------------------------------------

# Mac OS X 10.11, GCC 5.3
gcc -O3 -Wl,-twolevel_namespace -undefined error -flax-vector-conversions -framework Accelerate -I/usr/local/include -arch x86_64 -fopenmp "Estimate.c" -o "Estimate" -L/usr/local/lib -lgsl -lnlopt -lm

# Execute program
./Estimate
