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

    mixed_dummy_vars* eps = nu_s->new_eps_bins(0.0, a0, af, ms, N);
    
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

        
    mixed_dummy_vars* new_eps = nu_s->new_eps_bins(0.0, x_value, af, nu_s->get_ms(), y_values->get_num_bins());
    
    dummy_vars* old_eps = y_values->get_eps();
    // should this be dummy_vars* or mixed_dummy_vars* - does this matter?
    
    double T_cm = 1.0/x_value; // x_value = a0 ? T_cm = 1 / a0
    double results[6];
    
    for(int i = 0; i < y_values->get_num_bins(); i++){
        double eps_value = y_values->get_eps_value(i);
        
        // Get f values -> stored in results
        y_values->interpolate_extrapolate(eps_value, T_cm, results);
        for (int j = 0; j < 6; j++){
            // Set the f values in the freqs object to the new values that match the new eps
            y_values->set_f_value(i, j, results[j]);
        }
    }
    
    // What is the delete y_values for?
    /*
    delete y_values;
    y_values = new_y_values;
    */
    return af;
    
    /*int num_lin = new_eps->get_length() - new_eps->get_num_gl();
    
    freqs_ntT* new_y_values = new freqs_ntT(new_eps, y_values->get_sterile_density(), y_values->get_time(), y_values->get_Temp(), false);

    int old_tail_new_index =new_eps->get_length()-1;
    for( ; old_tail_new_index >= 0; old_tail_new_index--)
        if(new_eps->get_value(old_tail_new_index) < old_eps->get_max_linspace())
            break;
            
    double interpolated_values[6];
    for(int i = old_tail_new_index; i < new_eps->get_length(); i++){
        y_values->interpolated_f_values(new_eps->get_value(i), interpolated_values);
        for(int j = 0; j < 6; j++)
            new_y_values->set_f_value(i, j, interpolated_values[j]);
    }        
    
    dep_vars** new_nu_dist = new dep_vars*[6];
    dep_vars** old_nu_dist = new dep_vars*[6];
    for(int j = 0; j < 6; j++){
        new_nu_dist[j] = new dep_vars(y_values->get_num_bins());
        old_nu_dist[j] = new dep_vars(y_values->get_num_bins());
        new_y_values->get_neutrino_distribution(j, new_nu_dist[j]);
        y_values->get_neutrino_distribution(j, old_nu_dist[j]);
    }
    
    double old_cdf, new_cdf;
    
    y_values->interpolated_f_values(new_eps->get_value(old_tail_new_index), interpolated_values);
        
    for(int j = 0; j < 6; j++){
        old_cdf = old_eps->partial_integrate_pow_end(num_lin, old_nu_dist[j], 3);
        
        
        double blah = 0.5 * (old_eps->get_max_linspace() - new_eps->get_value(old_tail_new_index)) * (interpolated_values[j] * pow(new_eps->get_value(old_tail_new_index), 3)+ old_nu_dist[j]->get_value(num_lin-1) * pow(old_eps->get_max_linspace(),3));
        old_cdf += blah;
        new_cdf = new_eps->partial_integrate_pow_end(old_tail_new_index, new_nu_dist[j], 3);
        
        new_nu_dist[j]->multiply_by(old_cdf/new_cdf);
    }
    
    int bin_above;
    for(int i = old_tail_new_index-1; i > new_eps->get_num_gel(); i--){
        bin_above = old_eps->bin_below(new_eps->get_value(i)) + 1;
        y_values->interpolated_f_values(new_eps->get_value(i), interpolated_values);
        for(int j = 0; j < 6; j++){
            old_cdf = old_eps->partial_integrate_pow_end(bin_above, old_nu_dist[j], 3);
            old_cdf += 0.5 * (old_eps->get_value(bin_above) - new_eps->get_value(i)) * (interpolated_values[j] * pow(new_eps->get_value(i), 3));
            
            new_cdf = new_eps->partial_integrate_pow_end(i+1, new_nu_dist[j], 3);
            
            new_nu_dist[j]->set_value(i, (old_cdf - new_cdf) / pow(new_eps->get_value(i), 3) / new_eps->get_weight(i));            
        }
    }
    
    for(int i = 0; i <= new_eps->get_num_gel(); i++){
        y_values->interpolated_f_values(new_eps->get_value(i), interpolated_values);
        for(int j = 0; j < 6; j++)
            new_nu_dist[j]->set_value(i, interpolated_values[j]);
    }
        
    for(int j = 0; j < 6; j++){
        new_y_values->set_neutrino_distribution(j, new_nu_dist[j]);
        
        delete new_nu_dist[j];
        delete old_nu_dist[j];
    }    
    
    delete[] new_nu_dist;
    delete[] old_nu_dist;
    delete new_eps;
    
    delete y_values;
    y_values = new_y_values;
    
    return af;*/
    
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

