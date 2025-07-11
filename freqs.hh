#ifndef __FREQS_HH__
#define __FREQS_HH__

#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"

using std::ostream;

class freqs_ntT : public dep_vars{
    protected:
        int num_bins;
        
        mixed_dummy_vars* eps;
        
    public:
        freqs_ntT(mixed_dummy_vars*, double, double, double, bool);
        freqs_ntT(freqs_ntT*);
        ~freqs_ntT();
        
        mixed_dummy_vars* get_eps();
        void new_eps(mixed_dummy_vars*);
        int get_num_bins();
        
        double get_sterile_density();
        double get_time();
        double get_time_sec();
        double get_Temp();
        
        double get_eps_value(int);

        double get_f_value(int, int); // get_f_value(int bin, int type)
        void get_neutrino_distribution(int, dep_vars*);
                
        void set_sterile_density(double);
        void set_time(double);
        void set_Temp(double);
        
        void print_eps_file(ostream&);
        
        void set_f_value(int, int, double); // set_f_value(int bin, int type, double value)
        void set_neutrino_distribution(int, dep_vars*);
        
        void interpolated_f_values(double, double*);
        void interpolated_f_values(double, int, double*);
        
        void interpolate_extrapolate(double, double, double*);
        
        double interpolated_f_value(double);
        double interpolated_f_value(double, int);
        
        void neutrino_energy_and_pressure(double, double*, double*);
        double neutrino_energy(double);
};

double fifth_order_fit(double, double*, double*);
double interpolate_log_fifth(double, double*, double*);
double linear(double, double, double, double, double);
double interpolate_log_linear(double, double, double, double, double);
#endif