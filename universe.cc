#include <cmath>

#include "universe.hh"
#include "dummy_dep_vars.hh"
#include "constants.hh"

#include <iostream>

particle::particle(double m_ptcl)
{
    m = m_ptcl;
}

particle::~particle()
{ ;}

double particle::mass()
{   return m;  }

electron_positron::electron_positron() : particle(_electron_mass_) 
{ 
    x = new gl_dummy_vars(50);
    y = new dep_vars(50);    
}

electron_positron::~electron_positron()
{
    delete x;
    delete y;  
}

double electron_positron::energy_f(double u, double T){
    return (2 / (pow(_PI_,2))) * (pow(u,2) * pow(T,3) * sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1));
}

double electron_positron::pressure_f(double u, double T){
    return 2 / (3 * pow(_PI_,2)) * pow(u,4) * pow(T,5) / sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1);
}

double electron_positron::energy_temp_f(double u, double T){
    return (2 * T / (pow(_PI_,2))) * pow(u,2) * (pow(u,2) * pow(T,2) + pow(m,2)) * exp(sqrt(pow(u,2) * pow(T,2) + pow(m,2)) / T) / pow(exp(sqrt(pow(u,2) * pow(T,2) + pow(m,2)) / T) + 1,2);
}

double electron_positron::pressure_temp_f(double u, double T){
    return 2 / (3 * pow(_PI_,2)) * pow(u,4) * pow(T,3) * exp(sqrt(pow(u,2) * pow(T,2) + pow(m,2)) / T) / pow(exp(sqrt(pow(u,2) * pow(T,2) + pow(m,2)) / T) + 1,2);
}

double electron_positron::energy(double T)
{
    for(int i = 0; i < 50; i++)
        y->set_value(i, energy_f(x->get_value(i), T));
        
    return x->integrate(y);
}

double electron_positron::pressure(double T)
{
    for(int i = 0; i < 50; i++)
        y->set_value(i, pressure_f(x->get_value(i), T));
        
    return x->integrate(y);
}

double electron_positron::drho_dT(double T)
{
    for(int i = 0; i < 50; i++)
        y->set_value(i, energy_temp_f(x->get_value(i), T));
        
    return x->integrate(y);
}

double electron_positron::dP_dT(double T)
{
    for(int i = 0; i < 50; i++)
        y->set_value(i, pressure_temp_f(x->get_value(i), T));
        
    return x->integrate(y);
}

photon::photon() : particle(0) { ;}

double photon::energy(double T)
{   return 2. * _PI_ * _PI_ / 30. * pow(T, 4); }

double photon::pressure(double T)
{   return energy(T) / 3.; }

double photon::drho_dT(double T)
{   return energy(T) * 4. / T; }

double photon::dP_dT(double T)
{   return pressure(T) * 4. / T; }

thermal_neutrinos::thermal_neutrinos() : particle (0) { ;}

double thermal_neutrinos::energy(double T)
{   return 6. * 7. / 8. * _PI_ * _PI_ / 30. * pow(T, 4); }

double thermal_neutrinos::pressure(double T)
{   return energy(T) / 3.; }

double thermal_neutrinos::drho_dT(double T)
{   return energy(T) * 4. / T; }

double thermal_neutrinos::dP_dT(double T)
{   return pressure(T) * 4. / T; }

universe::universe(bool nu_th)
{
    N_particles = 2;
    
    if(nu_th)
        N_particles++;
        
    thermal_particles = new particle*[N_particles];
    thermal_particles[0] = new electron_positron;
    thermal_particles[1] = new photon;
    
    if(nu_th)
        thermal_particles[2] = new thermal_neutrinos;
}   

universe::~universe()
{
    for(int i = 0; i < N_particles; i++)
        delete thermal_particles[i];
    delete[] thermal_particles;
}

void universe::energy_and_pressure(double T, double* rho, double* P)
{
    *rho = 0.;
    *P = 0.;
    
    for(int i = 0; i < N_particles; i++){
        *rho += thermal_particles[i]->energy(T);
        *P += thermal_particles[i]->pressure(T);
    }
    return;
}

void universe::energy_pressure_and_derivs(double T, double* rho, double* P, double* drhodT, double* dPdT)
{
    *rho = 0.;
    *P = 0.;
    *drhodT = 0.;
    *dPdT = 0.;
    
    for(int i = 0; i < N_particles; i++){
        *rho += thermal_particles[i]->energy(T);
        *P += thermal_particles[i]->pressure(T);
        *drhodT += thermal_particles[i]->drho_dT(T);
        *dPdT += thermal_particles[i]->dP_dT(T);
    }
    return;
}

void universe::energy_entropy_and_derivs(double T, double* rho, double* s, double* dsdT){
    *rho = 0.;
    *s = 0.;
    *dsdT = 0.;
    
    double en = 0;
    double pr = 0;
    
    for(int i = 0; i < N_particles; i++){
        en = thermal_particles[i]->energy(T);
        pr = thermal_particles[i]->pressure(T);
        
        *rho += en;
        *s += (en + pr) / T;
        
        *dsdT += (thermal_particles[i]->drho_dT(T) + thermal_particles[i]->dP_dT(T)) / T;
        *dsdT += - (en+pr)/T / T;
    }
    return;   
}

