#include "collisions_MPI.hh"
#include "collisions.hh"
#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"
#include "freqs.hh"
#include "mpi.h"

#include <iostream>

using std::cout;
using std::endl;

collisions::collisions(int rank, int num_ranks, mixed_dummy_vars* e){
    myid = rank;
    numprocs = num_ranks;
    
    num_integrators = 0;
    
    eps = new mixed_dummy_vars(e);
    
    if(myid != 0){
        for(int i = myid-1; i < eps->get_length(); i += numprocs-1)
            num_integrators++;
        
        integrators = new collision_integral*[num_integrators];
        for(int i = myid-1; i < eps->get_length(); i += numprocs-1)
            integrators[i] = new nu_nu_collision(i, eps);
    }
}

collisions::~collisions(){
    delete eps;
    
    if(myid != 0){    
        for(int j = 0; j < num_integrators; j++)
            delete integrators[j];
        delete[] integrators;
    }
}

void collisions::compute_R(double Tcm, double T, freqs_ntT* output){
    output->zeros();
    
    double* out_vals = new double[6 * eps->get_length()]();
    double my_ans = 0.;
    int sender, tag;
    MPI_Status status;
    
    double dummy_int[6];
    
    if(myid == 0){
        MPI_Recv(dummy_int, 6, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        sender = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        
        for(int j = 0; j < 6; j++)
            out_vals[j*eps->get_length() + tag] += dummy_int[j];
        
        cout << "Received: " << tag << endl;
    }
    else{
        for(int j = 0; j < num_integrators; j++){
            integrators[j]->compute_R(Tcm, T, dummy_int);
            MPI_Send(dummy_int, 6, MPI_DOUBLE, 0, integrators[j]->get_bin(), MPI_COMM_WORLD);
        }
    }
    
    MPI_Bcast(out_vals, 6 * eps->get_length(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int i = 0; i < 6 * eps->get_length(); i++)
        output->set_value(i, out_vals[i]);
    
    delete[] out_vals;
}