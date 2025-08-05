#!/usr/bin/bash

rm coll

mpic++ run_collisions_MPI.cc collisions_MPI.cc collisions.cc mixed_dummy_vars.cc dummy_dep_vars.cc freqs.cc sterile.cc universe.cc -std=c++11 -o coll

mpiexec -n 6 coll