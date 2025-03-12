#include <cmath>

#include "freqs.hh"
#include "dummy_dep_vars.hh"
#include "constants.hh"

freqs_ntT::freqs_ntT(freqs_ntT* copy_me) : dep_vars(copy_me){
    num_bins = copy_me->get_num_bins();
    eps = new dummy_vars(copy_me->get_eps());
    
}

freqs_ntT::freqs_ntT(dummy_vars* e, double n0, double time, double Temp, bool thermal_dist) : dep_vars(6 * e->get_length() + 3){
    num_bins = e->get_length();
    eps = new dummy_vars(e);
    
    values[6 * num_bins] = n0;
    values[6 * num_bins + 1] = time;
    values[6 * num_bins + 2] = Temp;
    
    if(thermal_dist){
        dep_vars* f = new dep_vars(num_bins);
        for (int i = 0; i < num_bins; i++)
            f->set_value(i, 1./(exp(eps->get_value(i)) + 1));
            
        for(int j = 0; j < 6; j++)
            for(int i = 0; i < num_bins; i++)
                values[j * num_bins + i] = f->get_value(i);
                
        delete f;
    }
}

freqs_ntT::~freqs_ntT(){
    delete eps;
}

int freqs_ntT::get_num_bins()
{   return num_bins; }

dummy_vars* freqs_ntT::get_eps()
{   return eps; }

double freqs_ntT::get_sterile_density()
{   return values[6 * num_bins]; }

double freqs_ntT::get_time()
{   return values[6 * num_bins + 1]; }

double freqs_ntT::get_time_sec()
{   return values[6 * num_bins + 1] * _hbar_; }

double freqs_ntT::get_Temp()
{   return values[6 * num_bins + 2]; }

double freqs_ntT::get_f_value(int bin, int type)
{   return values[type * num_bins + bin];  }

void freqs_ntT::set_f_value(int bin, int type, double val)
{   values[type * num_bins + bin] = val; }

void freqs_ntT::set_sterile_density(double n0)
{   values[6 * num_bins] = n0; }

void freqs_ntT::set_time(double t0)
{   values[6 * num_bins + 1] = t0; }

void freqs_ntT::set_Temp(double T0)
{   values[6 * num_bins + 2] = T0; }

double freqs_ntT::neutrino_energy(double a)
{  
    double Tcm = 1./a;
    double coeff = pow(Tcm, 4) / (2 * _PI_ * _PI_);
    
    dep_vars* u = new dep_vars(num_bins);
    
    double rho = 0.;
    for(int j = 0; j < 6; j++){
        for(int i = 0; i < num_bins; i++)
            u->set_value(i, pow(eps->get_value(i), 3) * get_f_value(i, j));

        rho += eps->integrate(u);
    }
    rho *= coeff;    
    delete u;
    
    return rho;
}

void freqs_ntT::neutrino_energy_and_pressure(double a, double* rho, double* P)
{
    double Tcm = 1. / a;
    
    double coeff = pow(Tcm, 4) / (2 * _PI_ * _PI_);
    
    dep_vars* u = new dep_vars(num_bins);
    
    *rho = 0.;
    for(int j = 0; j < 6; j++){
        for(int i = 0; i < num_bins; i++)
            u->set_value(i, pow(eps->get_value(i), 3) * get_f_value(i, j));
        *rho += eps->integrate(u);
    }
    *rho *= coeff;
    *P = *rho / 3.;
    
    delete u;
}

void freqs_ntT::set_neutrino_distribution(int type, dep_vars* f){
    for(int i = 0; i < num_bins; i++)
        set_f_value(i, type, f->get_value(i));
}