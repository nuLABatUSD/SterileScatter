#include <cmath>
#include <iostream>
#include <fstream>

#include "sterile_decay.hh"
#include "sterile.hh"
#include "universe.hh"
#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "constants.hh"

using std::ostream;

sterile_decay::sterile_decay(double ms, double theta, double T0, mixed_dummy_vars* e) : ODESolve(){
    thermal = new universe(false);
    nu_s = new sterile(ms, theta);
    
    n0 = 0.174 * 3 * _zeta_3_ / 2. / _PI_ / _PI_ * pow(T0, 3);
    
    x_value = 1./T0;
    dx_value = 0.001 * x_value;
    
    y_values = new freqs_ntT(e, n0, 0., T0, true);
    
    std::cout << "m_s = " << ms << "MeV, lifetime (in seconds) = " << nu_s->get_lifetime_s() << std::endl;
}

sterile_decay::sterile_decay(double ms, double theta, double a0, double af, int N) : ODESolve(){
    thermal = new universe(false);
    nu_s = new sterile(ms, theta);

    mixed_dummy_vars* eps = nu_s->new_eps_bins(a0, af, ms, N); 
    // help! - should this call the original mdv constructor or the shift constructor?
    
    double T0 = 1./a0;

    n0 = 0.174 * 3 * _zeta_3_ / 2. / _PI_ / _PI_ * pow(T0, 3);
    
    x_value = 1./T0;
    dx_value = 0.001 * x_value;
    
    y_values = new freqs_ntT(eps, n0, 0., T0, true);
    
    std::cout << "m_s = " << ms << "MeV, lifetime (in seconds) = " << nu_s->get_lifetime_s() << std::endl;
    
    delete eps;
}

sterile_decay::~sterile_decay()
{   
    delete thermal;
    delete nu_s; }
    
/*double sterile_decay::shift_eps_by_multiple(double a_mult)
{   return shift_eps(x_value * a_mult);  }*/
    
double sterile_decay::shift_eps(double af)
{
    if(!nu_s->is_decay_on()){
        return af;}

    if(y_values->get_sterile_density() < n0 * fraction_n0_turn_off){
        std::cout << "Sterile Decays turned off" << std::endl;
        nu_s->turn_decay_off();
        return af;   
    }
    
    // ------------------------------------------------------
    
    freqs_ntT* new_f = new freqs_ntT(y_values);
    
    mixed_dummy_vars* new_eps = nu_s->new_eps_bins(y_values->get_a0(), x_value, af, nu_s->get_ms(), y_values->get_num_bins());
    // update to call the shift constructor
    
    new_f->new_eps(new_eps);
    
    double T_cm = 1.0/x_value;
    double results[6];
    int bin0;
    
    for(int i = 0; i < y_values->get_num_bins(); i++){
        double eps_value = new_eps->get_value(i);
        bin0 = y_values->interpolate_extrapolate(eps_value, T_cm, results);
        
        for (int j = 0; j < 6; j++){
            new_f->set_f_value(i, j, results[j]);
        }
        
        std::cout << eps_value << ", " << bin0 << ", " << y_values->get_eps_value(bin0) << std::endl;
    }
    
    delete y_values;
    y_values = new_f;
    return af;
    
}

bool sterile_decay::check_decay_on(){
    if(nu_s->is_decay_on() == true){
        return true;}
    else{
        return false;}
}
    
void sterile_decay::f(double a, freqs_ntT* inputs, freqs_ntT* derivs){
    dep_vars** p_all = new dep_vars*[6];
    for (int i = 0; i < 6; i++)
        p_all[i] = new dep_vars(inputs->get_num_bins());
        
/*    double rho_em, P_em, drhodT_em, dPdT_em;
    thermal->energy_pressure_and_derivs(inputs->get_Temp(), &rho_em, &P_em, &drhodT_em, &dPdT_em); */
    
    double rho_em, s_em, ds_em_dT;
    thermal->energy_entropy_and_derivs(inputs->get_Temp(), &rho_em, &s_em, &ds_em_dT);
        
    double ns = inputs->get_sterile_density();
    
    double rho = rho_em;
    rho += inputs->neutrino_energy(a);
    rho += nu_s->mass() * ns;
    
    double dtda = sqrt(3. / 8. / _PI_ / rho) * _planck_mass_ / a;
        
    nu_s->compute_full_term(inputs->get_eps(), 1./a, ns, dtda, p_all);
    
    double dQda = nu_s->get_rate() * nu_s->mass() * ns * a*a*a * dtda;

    dep_vars* eps_cubed = new dep_vars(inputs->get_num_bins());
    dep_vars* integrand = new dep_vars(inputs->get_num_bins());
    for(int j = 0; j < inputs->get_num_bins(); j++)
        eps_cubed->set_value(j, pow(inputs->get_eps_value(j), 3));
        
    double dq_neutrino = 0.;
    for(int i = 0; i < 6; i++){
        integrand->copy(p_all[i]);
        integrand->multiply_by(eps_cubed);
        dq_neutrino += inputs->get_eps()->integrate(integrand);
    }
    
    dq_neutrino *= (1./a) / (2. * _PI_ * _PI_);
    
    dQda -= dq_neutrino;
    
    double T = inputs->get_Temp();
    
    double dTda = dQda / T - 3 * a * a * s_em;
    dTda /= a * a * a;
    dTda /= ds_em_dT;
    
    double dnda = - nu_s->get_rate() * ns * dtda - 3 * ns / a;
    
    derivs->set_sterile_density(dnda);
    derivs->set_time(dtda);
    derivs->set_Temp(dTda);
    
    for(int i = 0; i < 6; i++)
        derivs->set_neutrino_distribution(i, p_all[i]);
    
    for(int i = 0; i < 6; i++)
        delete p_all[i];
    delete[] p_all;
    
    delete eps_cubed;
    delete integrand;
}

double sterile_decay::em_entropy_density(){
    double rho_em, s_em, ds_em_dT;
    thermal->energy_entropy_and_derivs(y_values->get_Temp(), &rho_em, &s_em, &ds_em_dT);
    
    return s_em;
}

void sterile_decay::print_eps_file(ostream& os)
{   y_values->print_eps_file(os);  }

double sterile_decay::get_Neff(){
    double Neff = 4./7. * pow(11./4., 4./3.) * 30. / (_PI_ * _PI_);
    Neff /= pow(y_values->get_Temp(), 4);
    Neff *= y_values->neutrino_energy(1./x_value);
    return Neff;
}   

freqs_ntT* sterile_decay::get_y_values(){
    return y_values;
}

double sterile_decay::get_neutrino_energy(double a){
    return y_values->neutrino_energy(a);
}

