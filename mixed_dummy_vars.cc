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


mixed_dummy_vars::mixed_dummy_vars(double a_0, double a_f, double m_s, int num) : dummy_vars(num){
// Update Protected Values
    a0 = a_0;
    af = a_f;
    ms = m_s;
    N = num;
    
// Comoving Temperatures
    double Tcm_0 = 1/a_0;
    double Tcm_f = 1/a_f;
    
// Find E-values
    calc_energies();
    
// Update max_linspace (inherited from dummy_vars)
    max_linspace = key_energies[5]/Tcm_f;
    
// Set Up Intervals
    intBounds[0] = 0.0;
    // Top Hat Intervals
    for(int i = 1; i < 5; i++){
        if(i % 2 == 1){ // if number is odd
            intBounds[i] = key_energies[i-1]/Tcm_0;}
        else{ // if number is even
            intBounds[i] = key_energies[i-1]/Tcm_f;}
    }
    // Delta Function Intervals
    int j = 0;
    for(int i = 5; i < 9; i++){
        if(i % 2 == 1){ // if number is odd
            intBounds[i] = key_energies[i-j-1]/Tcm_0;}
        else{ // if number is even
            intBounds[i] = key_energies[i-j-2]/Tcm_f;
            j++;} 
    }
    
    // Assign weights and values?
    int ind = 0;
    int gel_num, lin_num;
    
    for(int int_i = 0; int_i < def_N_ints; int_i++){
        // Set # of points per interval
        if(int_i == 0) gel_num = def_N_gel1;
        else gel_num = def_N_gel;
        
        if(int_i < 2) lin_num = def_N_linTop;
        else lin_num = def_N_linDel;
        
        // Create gel and linspace objects
        gel_dummy_vars* gel_vars = new gel_dummy_vars(gel_num, intBounds[int_i*2], intBounds[int_i*2+1]);        
        
        linspace_for_trap* lin_vars = new linspace_for_trap(intBounds[int_i*2+1], intBounds[int_i*2+2], lin_num);
    
        // for loops to set the weights and values from the _vars objects
        for(int i = 0; i < gel_num; i++){
            values[ind] = gel_vars->get_value(i);
            weights[ind] = gel_vars->get_weight(i);
            ind++;        
        }
        for(int i = 0; i < lin_num; i++){
            values[ind] = lin_vars->get_value(i);
            weights[ind] = lin_vars->get_weight(i);
            ind++; 
        }
        // delete the _vars 
        delete gel_vars;
        delete lin_vars;
    }
    
    int num_gl = N - ind;
    gl_dummy_vars* gl_vars = new gl_dummy_vars(num_gl, intBounds[8]);
    for(int i = 0; i < num_gl; i++){
        values[ind] = gl_vars->get_value(i);
        weights[ind] = gl_vars->get_weight(i);
        ind++;
    }
    
    delete gl_vars;  
}

mixed_dummy_vars::mixed_dummy_vars(double old_a0, double old_af, double new_af, double m_s, int num) : dummy_vars(num){
// Update Protected Values
    a0 = old_af;
    af = new_af;
    ms = m_s;
    N = num;
    
// Comoving Temperatures
    double old_Tcm_0 = 1/old_a0;
    double Tcm_0 = 1/a0;
    double Tcm_f = 1/af;
    
// Find E-values
    calc_energies();
    
// Update max_linspace (inherited from dummy_vars)
    max_linspace = key_energies[5]/Tcm_f;
    
// Set Up Intervals
    intBounds[0] = 0.0;
    // Top Hat Intervals
    for(int i = 1; i < 5; i++){
        if(i % 2 == 1){ // if number is odd
            intBounds[i] = key_energies[i-1]/Tcm_0;}
        else{ // if number is even
            intBounds[i] = key_energies[i-1]/Tcm_f;}
    }
    // Delta Function Intervals // UPDATE HERE - done?
    int j = 0;
    double key_energy;
    for(int i = 5; i < 11; i++){
        if(i < 8) key_energy = key_energies[4];
        else key_energy = key_energies[5];
        
        if(j == 0){ 
            intBounds[i] = key_energy/old_Tcm_0;
            j++;}
        else if (j == 1){
            intBounds[i] = key_energy/Tcm_0;
            j++;}
        else{ // if number is even
            intBounds[i] = key_energy/Tcm_f;
            j=0;} 
    }
    
    // Assign weights and values
    int ind = 0;
    int gel_num, lin_num1, lin_num2;
    
    for(int int_i = 0; int_i < def_N_ints; int_i++){
        // Set # of points per interval
        if(int_i < 2){
            if(int_i == 0) gel_num = def_N_gel1;
            else gel_num = def_N_gel;
            
            lin_num1 = def_N_linTop;
            
            // Create gel and linspace objects
            gel_dummy_vars* gel_vars = new gel_dummy_vars(gel_num, intBounds[int_i*2], intBounds[int_i*2+1]);        
            
            linspace_for_trap* lin_vars1 = new linspace_for_trap(intBounds[int_i*2+1], intBounds[int_i*2+2], lin_num1);
        
            // for loops to set the weights and values from the _vars objects
            for(int i = 0; i < gel_num; i++){
                values[ind] = gel_vars->get_value(i);
                weights[ind] = gel_vars->get_weight(i);
                ind++;        
            }
            for(int i = 0; i < lin_num1; i++){
                values[ind] = lin_vars1->get_value(i);
                weights[ind] = lin_vars1->get_weight(i);
                ind++; 
            }
            // delete the _vars 
            delete gel_vars;
            delete lin_vars1;
        }
        else{
            gel_num = def_N_gel;
            lin_num1 = 5; 
            lin_num2 = def_N_linDel - lin_num1;
            
            int bound_num;
            if(int_i == 2) bound_num = 4;
            else if(int_i == 3) bound_num = 7; 
            
            gel_dummy_vars* gel_vars = new gel_dummy_vars(gel_num, intBounds[bound_num], intBounds[bound_num+1]);        
            
            linspace_for_trap* lin_vars1 = new linspace_for_trap(intBounds[bound_num+1], intBounds[bound_num+2], lin_num1+1);
            linspace_for_trap* lin_vars2 = new linspace_for_trap(intBounds[bound_num+2], intBounds[bound_num+3], lin_num2);
            
            // for loops to set the weights and values from the _vars objects
            for(int i = 0; i < gel_num; i++){
                values[ind] = gel_vars->get_value(i);
                weights[ind] = gel_vars->get_weight(i);
                ind++;        
            }
            for(int i = 0; i < lin_num1; i++){
                values[ind] = lin_vars1->get_value(i);
                weights[ind] = lin_vars1->get_weight(i);
                ind++; 
            }
            for(int i = 0; i < lin_num2; i++){
                values[ind] = lin_vars2->get_value(i);
                weights[ind] = lin_vars2->get_weight(i);
                ind++;
            }
            // delete the _vars 
            delete gel_vars;
            delete lin_vars1;
            delete lin_vars2;
        }
    }
    
    int num_gl = N - ind;
    gl_dummy_vars* gl_vars = new gl_dummy_vars(num_gl, intBounds[10]);
    for(int i = 0; i < num_gl; i++){
        values[ind] = gl_vars->get_value(i);
        weights[ind] = gl_vars->get_weight(i);
        ind++;
    }
    delete gl_vars;
}

mixed_dummy_vars::mixed_dummy_vars(mixed_dummy_vars* copy_me) : dummy_vars(copy_me)
{
    // Update Protected Values?
    a0 = copy_me->get_a0();
    //double af;
    //double ms;
    //int N;
    
    for(int i = 0; i<6; i++){
        key_energies[i] = copy_me->get_key_energies(i);}
}

double mixed_dummy_vars::get_ms(){
    return ms;
}

double mixed_dummy_vars::get_a0(){
    return a0;
}

double mixed_dummy_vars::get_key_energies(int ind){
    return key_energies[ind];
}

double mixed_dummy_vars::get_intBounds(int ind){
    return intBounds[ind];
}

void mixed_dummy_vars::calc_energies(){
// Decays 1 and 2
    double E0_1 = ms / 2.0;
    double E0_2 = get_monoenergy(ms, 0.0, _neutral_pion_mass_);
    
// Decays 3 and 4
    double E_nu_prime = get_monoenergy(_charged_pion_mass_,0.0,_muon_mass_);
    double p_nu_prime = E_nu_prime;
    
    // Decay 3
    double gamma_3, v_3, p_3;
    compute_kinetics(ms, _charged_pion_mass_, _electron_mass_, &gamma_3, &v_3, &p_3);
    
    double Emin_3 = gamma_3 * (E_nu_prime - (v_3 * p_nu_prime));
    double Emax_3 = gamma_3 * (E_nu_prime + (v_3 * p_nu_prime));
    
    // Decay 4
    double gamma_4, v_4, p_4;
    compute_kinetics(ms, _charged_pion_mass_, _muon_mass_, &gamma_4, &v_4, &p_4);
    
    double Emin_4 = gamma_4 * (E_nu_prime - (v_4 * p_nu_prime));
    double Emax_4 = gamma_4 * (E_nu_prime + (v_4 * p_nu_prime));
    
    // Update energies array (protected in mixed_dummy_vars)
    key_energies[0] = Emin_3;
    key_energies[1] = Emin_4;
    key_energies[2] = Emax_4;
    key_energies[3] = Emax_3;
    key_energies[4] = E0_2;
    key_energies[5] = E0_1;

}
