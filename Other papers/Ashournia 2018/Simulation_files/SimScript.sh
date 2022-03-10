#!/bin/sh
#----------------------------------------------------
# By Damoun Ashournia, University of Oxford
#----------------------------------------------------
# Bash script to compile and run the simulation 
# module.
#----------------------------------------------------

# Mac OS X 10.10, GCC 4.9
gcc -O3 -Wl,-twolevel_namespace -undefined error -flax-vector-conversions -framework Accelerate -Wall -I/usr/local/include -arch x86_64 -fopenmp "Simulate.c" -o "Simulate" -L/usr/local/lib -lgsl

# Execute program
./Simulate
