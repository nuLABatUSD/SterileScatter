#ifndef _COLLISIONS_MPI_HH_
#define _COLLISIONS_MPI_HH_

#include "collisions.hh"
#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"
#include "freqs.hh"

class collisions{
    protected:
        int myid, numprocs, num_integrators;
        
        mixed_dummy_vars* eps;
        collision_integral** integrators;
        
    public:
        collisions(int, int, mixed_dummy_vars*);
        ~collisions();
        
        void compute_R(double, double, freqs_ntT*);
        void collision_integral(freqs_ntT*, bool, freqs_ntT*);
}


#endif