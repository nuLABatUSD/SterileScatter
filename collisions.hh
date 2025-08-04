#ifndef _COLLISIONS_HH_
#define _COLLISIONS_HH_

#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "mixed_dummy_vars.hh"

#include <iostream>

using std::cout;
using std::endl;

const double _COMPUTE_R_ERROR_ = -999.;

class collision_integral{
    protected:
        int bin, N_bins;
        double eps_value;
        
        mixed_dummy_vars* eps;
        
        dep_vars* outer_vals;
        dummy_vars* outer_dummy_vars;
        
        dep_vars** inner_vals;
        dummy_vars** inner_dummy_vars;
        
        double*** F_values;

    public:
        collision_integral(int, mixed_dummy_vars*); //(int bin, int flavor, dummy_vars* eps);
        virtual ~collision_integral() = default;
        
        int get_bin();
        int get_num_bins();
        double get_eps_value();
        
        virtual void populate_F(freqs_ntT*, double, bool) = 0;
        virtual double interior_integral(double, int, int) = 0;
        virtual void whole_integral(freqs_ntT*, double, bool, double*) = 0;
        
        virtual void compute_R(double, double, double*) = 0;
};

class nu_nu_collision : public collision_integral{

    public:
        nu_nu_collision(int, mixed_dummy_vars*);
        ~nu_nu_collision();
        
        void populate_F(freqs_ntT*, double, bool);
        double interior_integral(double, int, int);
        void whole_integral(freqs_ntT*, double, bool, double*);
        
        void compute_R(double, double, double*);
        
        double J(double, double, double);


};

double J1(double, double, double);
double J2(double, double);
double J3(double, double, double);

#endif