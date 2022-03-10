#!/bin/sh
#--------------------------------------------------------
# By Damoun Ashournia, University of Oxford
#--------------------------------------------------------
# Bash script to compile and run the std. err. module.
#--------------------------------------------------------

# Mac OS X 10.10, GCC4.9
gcc -O3 -Wl,-twolevel_namespace -undefined error -flax-vector-conversions -framework Accelerate -I/usr/local/include -arch x86_64 -fopenmp "StdErrors.c" -o "StdErrors" -L/usr/local/lib -lgsl -lm

# Execute program
./StdErrors