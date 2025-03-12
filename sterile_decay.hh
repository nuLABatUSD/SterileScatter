#ifndef __STERILE_DECAY_HH__
#define __STERILE_DECAY_HH__

#include "universe.hh"
#include "sterile.hh"
#include "freqs.hh"
#include "ODESolve.hh"

class sterile_decay : public ODESolve<freqs_ntT>{
    protected:
        universe* thermal;
        sterile* nu_s;
        dummy_vars* eps;
        
    public:
        sterile_decay(double, double, double, dummy_vars*);
        ~sterile_decay();

        void f(double, freqs_ntT*, freqs_ntT*);
        
        double em_entropy_density();
};

#endif