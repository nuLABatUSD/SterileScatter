#ifndef __STERILE_HH__
#define __STERILE_HH__

#include "universe.hh"

const double sigma_delta_boxcar = 0.05;

class sterile : public particle{
    protected:
        double sin2_2th;
        
        double decay_rate;
        double lifetime_MeV;
        double lifetime_s;
        
        double rate[5];
        
    public:
        sterile(double, double);
        
        double get_rate();
        double get_lifetime();
        double get_lifetime_s();
        
        void calc_rates();
        
        double get_decay_type_one(double, double, double);
        double get_decay_type_two(double, double);
        double get_decay_type_three(double);
        double get_decay_type_four(double, double);
        
        void compute_dPdtdE(dummy_vars*, double, dep_vars**);
        void compute_full_term(dummy_vars*, double, dep_vars**);
        
        double min_low();
        
        double energy(double);
        double pressure(double);
        double drho_dT(double);
        double dP_dT(double);
       
};

void compute_kinetics(double, double, double, double*, double*, double*);
double get_monoenergy(double, double, double);

double get_gamma_a(double);
double get_gamma_b(double);
double type_four_integrand(double, double);


#endif