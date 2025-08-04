#include "collisions.hh"
#include "mixed_dummy_vars.hh"
#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "collisions_MPI.hh"
#include "mpi.h"

#include <iostream>
#include <chrono>

using std::cout;
using std::endl;

using namespace std::chrono;

int main(int argc, char* argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    mixed_dummy_vars* eps = new mixed_dummy_vars(0.10, 0.11, 300., 200);

    collisions* C_MPI = new collisions(myid, numprocs, eps);
    
    freqs_ntT* R_values = new freqs_ntT(eps, 0., 0., 10., false);

    auto start = high_resolution_clock::now();
    
    C_MPI->compute_R(10., 10., R_values);
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    if(myid == 0)
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;

    delete R_values;
    delete C_MPI;
        
    delete eps;
    
    MPI_Finalize();   
}