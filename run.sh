#!/usr/bin/bash

rm test
#rm outputs/sim_output*
#rm outputs/eps_output*
#rm outputs/ned_output*

g++ sterile.cc sterile_decay.cc dummy_dep_vars.cc mixed_dummy_vars.cc freqs.cc universe.cc run_main.cc -std=c++11 -o test

./test 300 5e-5 test_1
./test 300 4e-5 test_2
./test 300 3e-5 test_3
./test 300 2e-5 test_4
./test 300 1e-5 test_5
