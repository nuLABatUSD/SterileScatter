#ifndef __FREQS_HH__
#define __FREQS_HH__

#include "dummy_dep_vars.hh"

class freqs_ntT : public dep_vars{
    protected:
        int num_bins;
        
        dummy_vars* eps;
        
    public:
        freqs_ntT(dummy_vars*, double, double, double, bool);
        freqs_ntT(freqs_ntT*);
        ~freqs_ntT();
        
        dummy_vars* get_eps();
        int get_num_bins();
        
        double get_sterile_density();
        double get_time();
        double get_time_sec();
        double get_Temp();
        
        double get_f_value(int, int);
        
        void set_sterile_density(double);
        void set_time(double);
        void set_Temp(double);
        
        void set_f_value(int, int, double);
        void set_neutrino_distribution(int, dep_vars*);
        
        void neutrino_energy_and_pressure(double, double*, double*);
        double neutrino_energy(double);
};

#endif