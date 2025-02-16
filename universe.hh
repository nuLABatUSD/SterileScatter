#ifndef __UNIVERSE_HH__
#define __UNIVERSE_HH__

#include "dummy_dep_vars.hh"

class particle;

class universe{
    protected:
        int N_particles;
        particle** thermal_particles;
    
    public:
        universe(bool);
        ~universe();    
        
        void energy_and_pressure(double, double*, double*);
        void energy_pressure_and_derivs(double, double*, double*, double*, double*);
};

class particle{
    protected:
        double m;
        
    public:
        particle(double);
        virtual ~particle();
        
        double mass();
        
        virtual double energy(double) = 0;
        virtual double pressure(double) = 0;
        virtual double drho_dT(double) = 0;
        virtual double dP_dT(double) = 0;
};

class electron_positron : public particle{
    protected:
        gl_dummy_vars* x;
        dep_vars* y;

    public:
        electron_positron();
        ~electron_positron();
        
        double energy_f(double, double);
        double pressure_f(double, double);
        
        double energy_temp_f(double, double);
        double pressure_temp_f(double, double);
        
        double energy(double);
        double pressure(double);
        double drho_dT(double);
        double dP_dT(double);
};

class photon : public particle{
    public:
        photon();
        
        double energy(double);
        double pressure(double);
        double drho_dT(double);
        double dP_dT(double);        
};

class thermal_neutrinos : public particle{
    public:
        thermal_neutrinos();
        
        double energy(double);
        double pressure(double);
        double drho_dT(double);
        double dP_dT(double);        
};

#endif