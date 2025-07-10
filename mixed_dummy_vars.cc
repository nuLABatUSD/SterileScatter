#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "mixed_dummy_vars.hh"
#include "constants.hh"
#include "sterile.hh"

#include "dummy_dep_vars.hh"
#include "gl_vals.hh"
#include "gel_vals.hh"

using std::cout;
using std::endl;
using std::ostream;


mixed_dummy_vars::mixed_dummy_vars(double E_start, double a_0, double a_f, double m_s, int num) : dummy_vars(num){
// Update Protected Values
    Estart = E_start;
    a0 = a_0;
    af = a_f;
    ms = m_s;
    N = num;
    
// Comoving Temperatures
    double Tcm_0 = 1/a_0;
    double Tcm_f = 1/a_f;
    
// Find E-values
    double E0_1 = m_s / 2.0;
    double E0_2 = get_monoenergy(m_s, 0.0, _neutral_pion_mass_);
    
    // Decays 3 and 4
    double E_nu_prime = get_monoenergy(_charged_pion_mass_,0.0,_muon_mass_);
    double p_nu_prime = E_nu_prime;
    
    // Decay 3
    double gamma_3, v_3, p_3;
    compute_kinetics(m_s, _charged_pion_mass_, _electron_mass_, &gamma_3, &v_3, &p_3);
    
    double Emin_3 = gamma_3 * (E_nu_prime - (v_3 * p_nu_prime));
    double Emax_3 = gamma_3 * (E_nu_prime + (v_3 * p_nu_prime));
    
    // Decay 4
    double gamma_4, v_4, p_4;
    compute_kinetics(m_s, _charged_pion_mass_, _muon_mass_, &gamma_4, &v_4, &p_4);
    
    double Emin_4 = gamma_4 * (E_nu_prime - (v_4 * p_nu_prime));
    double Emax_4 = gamma_4 * (E_nu_prime + (v_4 * p_nu_prime));
    
// Update max_linspace (inherited from dummy_vars)
    max_linspace = E0_1/Tcm_f;
    
// Make first 10 Gauss-Legendre points
    gel_dummy_vars* gel_1 = new gel_dummy_vars(10, E_start/Tcm_0, Emin_3/Tcm_0);
    
// Make Linspace points for Emin_3 and Emin_4
    linspace_for_trap* lin_1 = new linspace_for_trap(Emin_3/Tcm_0, Emin_4/Tcm_f, 50);

// Make next Gauss-Legendre points
    gel_dummy_vars* gel_2 = new gel_dummy_vars(5, Emin_4/Tcm_f, Emax_4/Tcm_0);
    
// Make Linspace  points for Emax_3 and Emax_4
    linspace_for_trap* lin_2 = new linspace_for_trap(Emax_4/Tcm_0, Emax_3/Tcm_f, 50);
    
// Make next Gauss-Legendre points
    gel_dummy_vars* gel_3 = new gel_dummy_vars(5, Emax_3/Tcm_f, E0_2/Tcm_0);
    
// Make Linspace points for E0_2
    linspace_for_trap* lin_3 = new linspace_for_trap(E0_2/Tcm_0, E0_2/Tcm_f, 35);
    
// Make next Gauss-Legendre points
    gel_dummy_vars* gel_4 = new gel_dummy_vars(5, E0_2/Tcm_f, E0_1/Tcm_0);
    
// Make Linspace points for E0_1
    linspace_for_trap* lin_4 = new linspace_for_trap(E0_1/Tcm_0, E0_1/Tcm_f, 35);
    
// Make 5 Gauss-Laguerre points
    gl_dummy_vars* gl = new gl_dummy_vars(5, E0_1/Tcm_f);

// Combine values and weights from each segment
    for(int i = 0; i < 10; i++){
        values[i] = gel_1->get_value(i);
        weights[i] = gel_1->get_weight(i);
    }
    for(int i = 10; i < 60; i++){
        values[i] = lin_1->get_value(i-10);
        weights[i] = lin_1->get_weight(i-10);
    }
    for(int i = 60; i < 65; i++){
        values[i] = gel_2->get_value(i-60);
        weights[i] = gel_2->get_weight(i-60);
    }
    for(int i = 65; i < 115; i++){
        values[i] = lin_2->get_value(i-65);
        weights[i] = lin_2->get_weight(i-65);
    }
    for(int i = 115; i < 120; i++){
        values[i] = gel_3->get_value(i-115);
        weights[i] = gel_3->get_weight(i-115);
    }
    for(int i = 120; i < 155; i++){
        values[i] = lin_3->get_value(i-120);
        weights[i] = lin_3->get_weight(i-120);
    }
    for(int i = 155; i < 160; i++){
        values[i] = gel_4->get_value(i-155);
        weights[i] = gel_4->get_weight(i-155);
    }
    for(int i = 160; i < 195; i++){
        values[i] = lin_4->get_value(i-160);
        weights[i] = lin_4->get_weight(i-160);
    }
    for(int i = 195; i < 200; i++){
        values[i] = gl->get_value(i-195);
        weights[i] = gl->get_weight(i-195);
    }
    
    delete gel_1;
    delete lin_1;
    delete gel_2;
    delete lin_2;
    delete gel_3;
    delete lin_3;
    delete gel_4;
    delete lin_4;
    delete gl;
    
    // Update energies array (protected in mixed_dummy_vars)
    key_energies[0] = Emin_3;
    key_energies[1] = Emin_4;
    key_energies[2] = Emax_4;
    key_energies[3] = Emax_3;
    key_energies[4] = E0_2;
    key_energies[5] = E0_1;
    
}

mixed_dummy_vars::mixed_dummy_vars(mixed_dummy_vars* copy_me)
{
    N = copy_me->get_length();
    values = new double[N]();
    weights = new double[N]();
    max_linspace = copy_me->get_max_linspace();

    for(int i = 0; i<N; i++)
        {
            values[i] = copy_me->get_value(i);
            weights[i] = copy_me->get_weight(i);
        }
}

double mixed_dummy_vars::get_key_energies(int ind){
    return key_energies[ind];
}