#ifndef __STERILE_HH__
#define __STERILE_HH__

#include "universe.hh"
#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"

const double sigma_delta_boxcar = 0.05;

class sterile : public particle{
    protected:
        double theta;
        double sin2_2th;
        
        double decay_rate;
        double lifetime_MeV;
        double lifetime_s;
        
        double rate[5];
        
        double E_low, E_high;
        double energies[6];
        bool decay_on;
        
        double ms;// Added by Albert
    public:
        sterile(double, double);
        sterile(sterile*);
        
        double get_ms(); // Added by Albert
        
        double get_theta();
        
        double get_rate();
        double get_lifetime();
        double get_lifetime_s();
        
        bool is_decay_on();
        void turn_decay_off();
        
        void calc_rates();
        void calculate_min_max_energy();
        void calculate_energies();
        
        double get_E_low();
        double get_E_high();
        
        gel_linspace_gl* new_eps_bins(double, double, int);
        mixed_dummy_vars* new_eps_bins(double, double, double, int);
        mixed_dummy_vars* new_eps_bins(double, double, double, double, int); // shift constructor
        
        double get_decay_type_one(double, double, double);
        double get_decay_type_two(double, double);
        double get_decay_type_three(double);
        double get_decay_type_four(double, double);
        
        void compute_dPdtdE(dummy_vars*, double, dep_vars**);
        void compute_full_term(dummy_vars*, double, double, double, dep_vars**);
        
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