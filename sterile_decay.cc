#include <cmath>
#include <iostream>

#include "sterile_decay.hh"
#include "sterile.hh"
#include "universe.hh"
#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "constants.hh"

sterile_decay::sterile_decay(double ms, double theta, double T0, dummy_vars* e) : ODESolve(){
    thermal = new universe(false);
    nu_s = new sterile(ms, theta);
    eps = new dummy_vars(e);
    
    double n0 = 0.174 * 3 * _zeta_3_ / 2. / _PI_ / _PI_ * pow(T0, 3);
    
    x_value = 1./T0;
    dx_value = 0.001 * x_value;
    
    y_values = new freqs_ntT(eps, n0, 0., T0, true);
    
    std::cout << "m_s = " << ms << "MeV, lifetime (in seconds) = " << nu_s->get_lifetime_s() << std::endl;
}

sterile_decay::sterile_decay(double ms, double theta, double a0, double af, int N) : ODESolve(){
    thermal = new universe(false);
    nu_s = new sterile(ms, theta);

    std::cout << "E_low = " << nu_s->get_E_low() << ", E_high = " << nu_s->get_E_high() << std::endl;
    eps = nu_s->new_eps_bins(a0, af, N);
    
    double T0 = 1./a0;

    double n0 = 0.174 * 3 * _zeta_3_ / 2. / _PI_ / _PI_ * pow(T0, 3);
    
    x_value = 1./T0;
    dx_value = 0.001 * x_value;
    
    y_values = new freqs_ntT(eps, n0, 0., T0, true);
    
    std::cout << "m_s = " << ms << "MeV, lifetime (in seconds) = " << nu_s->get_lifetime_s() << std::endl;
    
    eps->print_all();
}

sterile_decay::~sterile_decay()
{   delete eps;
    delete thermal;
    delete nu_s; }
    
void sterile_decay::f(double a, freqs_ntT* inputs, freqs_ntT* derivs){
    dep_vars** p_all = new dep_vars*[6];
    for (int i = 0; i < 6; i++)
        p_all[i] = new dep_vars(eps->get_length());
        
/*    double rho_em, P_em, drhodT_em, dPdT_em;
    thermal->energy_pressure_and_derivs(inputs->get_Temp(), &rho_em, &P_em, &drhodT_em, &dPdT_em); */
    
    double rho_em, s_em, ds_em_dT;
    thermal->energy_entropy_and_derivs(inputs->get_Temp(), &rho_em, &s_em, &ds_em_dT);
        
    double ns = inputs->get_sterile_density();
    
    double rho = rho_em;
    rho += inputs->neutrino_energy(a);
    rho += nu_s->mass() * ns;
    
    double dtda = sqrt(3. / 8. / _PI_ / rho) * _planck_mass_ / a;
        
    nu_s->compute_full_term(eps, 1./a, ns, dtda, p_all);
    
    double dQda = nu_s->get_rate() * nu_s->mass() * ns * a*a*a * dtda;

    dep_vars* eps_cubed = new dep_vars(eps->get_length());
    dep_vars* integrand = new dep_vars(eps->get_length());
    for(int j = 0; j < eps->get_length(); j++)
        eps_cubed->set_value(j, pow(eps->get_value(j), 3));
        
    double dq_neutrino = 0.;
    for(int i = 0; i < 6; i++){
        integrand->copy(p_all[i]);
        integrand->multiply_by(eps_cubed);
        dq_neutrino += eps->integrate(integrand);
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
}

double sterile_decay::em_entropy_density(){
    double rho_em, s_em, ds_em_dT;
    thermal->energy_entropy_and_derivs(y_values->get_Temp(), &rho_em, &s_em, &ds_em_dT);
    
    return s_em;
}