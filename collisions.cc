#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "mixed_dummy_vars.hh"
#include "collisions.hh"
#include "constants.hh"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

using std::abs;

collision_integral::collision_integral(int b, mixed_dummy_vars* e){
    bin = b;
    N_bins = e->get_length();
    eps_value = e->get_value(bin);
    
    eps = new mixed_dummy_vars(e);
}

collision_integral::~collision_integral(){
    delete eps;
}

int collision_integral::get_bin()
{   return bin;}

int collision_integral::get_num_bins()
{   return N_bins;}

double collision_integral::get_eps_value()
{   return eps_value;}

nu_nu_collision::nu_nu_collision(int b, mixed_dummy_vars* e) : collision_integral(b, e){
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

void nu_nu_collision::populate_F(freqs_ntT* f, double Tcm, bool net){
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
            p3_value = inner_dummy_vars[p2]->get_value(p3);
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

double J1(double p1, double p2, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1+p2,2) - 15 * (p1+p2) * p3 + 6 * pow(p3,2));  
}

double J2(double p1, double p2){
    return 16./15 * pow(p2,3) * (10 * pow(p1,2) + 5 * p1*p2 + pow(p2,2));  
}

double J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2,5) - 10 * pow(p1+p2,2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3,4) - 6 * pow(p3,5));
}

double nu_nu_collision::J(double p1, double p2, double p3){
    if(p2 < p1){
        if(p3 < p2)
            return J1(p1, p2, p3);
        else if (p3 < p1)
            return J2(p1, p2);
        else if (p3 < p1 + p2)
            return J3(p1, p2, p3);
        else
            return 0;
    }
    else{
        if(p3 < p1)
            return J1(p1, p2, p3);
        else if (p3 < p2)
            return J2(p2, p1);
        else if (p3 < p1 + p2)
            return J3(p1, p2, p3);
        else
            return 0;
    }
}

double nu_nu_collision::interior_integral(double Tcm, int p2, int which_term){
    inner_vals[p2]->zeros();
    
    for(int p3 = 0; p3 < inner_vals[p2]->get_length(); p3++)
        inner_vals[p2]->set_value(p2, J1(eps_value, outer_dummy_vars->get_value(p2), inner_dummy_vars[p2]->get_value(p3)) * F_values[which_term][p2][p3]);
    
    return inner_dummy_vars[p2]->integrate(inner_vals[p2]);
}

void nu_nu_collision::whole_integral(freqs_ntT* f, double Tcm, bool net, double* results){
    if (eps_value == 0)
        for(int j = 0; j < 6; j++)
            results[j] = 0.;
    else{
        populate_F(f, Tcm, net);
        
        double coeff = pow(Tcm, 5) * pow(_GF_, 2) / (pow(2*_PI_, 3) * pow(eps_value,2 ));
        
        for(int j = 0; j < 6; j++){
            for(int p2 = 0; p2 < N_bins; p2++)
                outer_vals->set_value(p2, interior_integral(Tcm, p2, j));
            results[j] = coeff * outer_dummy_vars->integrate(outer_vals);
        }
    }
}

void nu_nu_collision::compute_R(double Tcm, double T, double* results){
    freqs_ntT* thermal = new freqs_ntT(eps, 0., 0., Tcm, true);
    
    double net[6];
    double frs[6];
    
    whole_integral(thermal, Tcm, true, net);
    whole_integral(thermal, Tcm, false, frs);
    
    for(int j = 0; j < 6; j++){
        if (net[j] + frs[j] == 0)
            results[j] = _COMPUTE_R_ERROR_;
        else
            results[j] = abs(net[j] / frs[j]);
       }
    
    delete thermal;
}


