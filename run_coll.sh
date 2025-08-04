#!/usr/bin/bash

rm test

g++ run_collisions.cc collisions.cc mixed_dummy_vars.cc dummy_dep_vars.cc freqs.cc sterile.cc universe.cc -std=c++11 -o test

./test