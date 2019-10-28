#!/usr/bin/env bash
set -e
gcc -Wall -O3 -I/home/juvid/gsl-2.5/include -c smc.c
gcc -L/home/juvid/gsl/lib smc.o -lgsl -lgslcblas -lm -o smc.ce
./smc.ce
