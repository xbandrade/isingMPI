#!/bin/bash

make clean
make
echo ""
mpirun -np 4 ./ising.x
