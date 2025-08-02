#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "mixed_dummy_vars.hh"
#include "collisions.hh"

#include <iostream>

using std::cout;
using std::endl;

collision_integral::collision_integral(int b, dummy_vars* e){
    bin = b;
    N_bins = e->get_length();
    eps_value = e->get_value(bin);
}

nu_nu_collision::nu_nu_collision(int b, dummy_vars* e) : collision_integral(b, e){
    outer_dummy_vars = new dummy_vars(e);
    outer_vals = new dep_vars(N_bins);
    
    inner_dummy_vars = new dummy_vars*[N_bins];
    inner_vals = new dep_vars*[N_bins];

    F_values = new double**[6];
    for(int j = 0; j < 6; j++)
        F_values[j] = new double*[N_bins];

    for(int i = 0; i < N_bins; i++){
        inner_dummy_vars[i] = new dummy_vars(e);
        inner_vals[i] = new dep_vars(N_bins);
        for(int j = 0; j < 6; j++)
            F_values[j][i] = new double[N_bins]();
    }
    
    
}

nu_nu_collision::~nu_nu_collision(){
    for(int i = 0; i < N_bins; i++){
        delete inner_dummy_vars[i];
        delete inner_vals[i];
        for(int j = 0; j < 6; j++)
            delete[] F_values[j][i];
    }
    
    for(int j = 0; j < 6; j++)
        delete[] F_values[j];
        
    delete[] F_values;
    
    delete[] inner_dummy_vars;
    delete[] inner_vals;

    delete outer_dummy_vars;
    delete outer_vals;
}

void nu_nu_collision::populate_F(freqs* f, double Tcm, bool net){
    double p2_value, p3_value, p3_max, p4_value;
    double interp_values[6];
    int target;
    
    double rev;
    if(net)
        rev = -1.;
    else
        rev = 1.;
        
    for(int p2 = 0; p2 < N_bins; p2++){
        p2_value = outer_dummy_vars->get_value(p2);
        p3_max = eps_value + p2_value;
        for(int p3 = 0; p3 < N_bins; p3++){
            p3_value = inner_dummy_vars->get_value(p3);
            if(p3_value > p3_max){
                for(int j = 0; j < 6; j++)
                    F_values[j][p2][p3] = 0.; }
            else{
                p4_value = p3_max - p3_value;
                
                f->interpolate_extrapolate(p4_value, Tcm, interp_values);
                
                for(int j = 0; j < 6; j++){
                    F_values[j][p2][p3] = rev * f->get_f_value(bin, j) * f->get_f_value(p2, j) * (1. - f->get_f_value(p3, j)) * (1. - interp_values[j]);
                    F_values[j][p2][p3] += (1. - f->get_f_value(bin, j)) * (1. - f->get_f_value(p2, j)) * f->get_f_value(p3, j) * interp_values[j];
                    
                    for(int k = 2; k < 5; k += 2){
                        target = (j + k) % 6;
                        F_values[j][p2][p3] += rev * 0.5 * f->get_f_value(bin, j) * f->get_f_value(p2, target) * (1. - f->get_f_value(p3, target)) * (1. - interp_values[j]);
                        F_values[j][p2][p3] += 0.5 * (1. - f->get_f_value(bin, j)) * (1. - f->get_f_value(p2, target)) * f->get_f_value(p3, target) * interp_values[j];
                    }
                }
            }
        }
    }
}

double interior_integral(double Tcm, int p2, int which_term){
    inner_vals[p2]->zeros();
    
    
    return inner_dummy_vars[p2]->integrate(inner_vals[p2]);
}



