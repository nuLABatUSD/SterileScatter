#!/usr/bin/bash

rm test

mpic++ run_collisions_MPI.cc collisions_MPI.cc collisions.cc mixed_dummy_vars.cc dummy_dep_vars.cc freqs.cc sterile.cc universe.cc -std=c++11 -o test

mpiexec -n 4 test