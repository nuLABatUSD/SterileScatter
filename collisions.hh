#ifndef _COLLISIONS_HH_
#define _COLLISIONS_HH_

#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "mixed_dummy_vars.hh"

#include <iostream>

using std::cout;
using std::endl;

class collision_integral{
    protected:
        int bin, N_bins;
        double eps_value;
        
        dep_vars* outer_vals;
        dummy_vars* outer_dummy_vars;
        
        dep_vars** inner_vals;
        dummy_vars** inner_dummy_vars;
        
        double*** F_values;

    public:
        collision_integral(int, int, dummy_vars*); //(int bin, int flavor, dummy_vars* eps);
        
        virtual void populate_F(freqs*, bool);
        virtual double interior_integral(int)
        virtual void whole_integral(freqs*, bool, double*);
};

class nu_nu_collision : public collision_integral{

    public:
        nu_nu_collision(int, int, dummy_vars*);
        
        void populate_F(freqs*, double, bool);
        double interior_integral(double, int, int);
        void whole_integral(freqs*, double, bool, double*);
        
        double J(double, double, double);


};

double J1(double, double, double);
double J2(double, double);
double J3(double, double, double);

#endif