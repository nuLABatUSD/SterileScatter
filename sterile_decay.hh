#ifndef __STERILE_DECAY_HH__
#define __STERILE_DECAY_HH__

#include "universe.hh"
#include "sterile.hh"
#include "freqs.hh"
#include "ODESolve.hh"

using std::ostream;

const double fraction_n0_turn_off = 1.e-6;

class sterile_decay : public ODESolve<freqs_ntT>{
    protected:
        universe* thermal;
        sterile* nu_s;
        
        double n0;
        
    public:
        sterile_decay(double, double, double, dummy_vars*);
        sterile_decay(double, double, double, double, int);
        ~sterile_decay();

        double shift_eps(double);
        double shift_eps_by_multiple(double);
        
        void f(double, freqs_ntT*, freqs_ntT*);
        
        double em_entropy_density();
        
        void print_eps_file(ostream&);
        
        double get_Neff();
};

#endif