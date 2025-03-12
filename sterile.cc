#include <cmath>
#include <iostream>

#include "sterile.hh"
#include "constants.hh"

sterile::sterile(double m_s, double theta) : particle(m_s){
    sin2_2th = pow(sin(theta), 2);
    
    calc_rates();
    
//    rate[1] = 0;
//    rate[2] = 0;
//    rate[3] = 0;
//    rate[4] = 0;
    
    decay_rate  = 3 * rate[1];
    decay_rate += 3 * rate[2];
    decay_rate += 2 * rate[3];
    decay_rate += 2 * rate[4];
        
    lifetime_MeV = 1./decay_rate;
    lifetime_s = lifetime_MeV * _hbar_;
}

double sterile::get_rate()
{   return decay_rate;  }

double sterile::get_lifetime()
{   return lifetime_MeV;    }

double sterile::get_lifetime_s()
{   return lifetime_s;  }

void sterile::calc_rates(){
    rate[1] = 3 * _fine_structure_ * pow(_GF_,2) * pow(m,5) * sin2_2th / (512 * pow(_PI_,4));
    
    if (m >= _neutral_pion_mass_)
        rate[2] = pow(_GF_,2) * pow(_pion_decay_,2) * m * (pow(m,2) - pow(_neutral_pion_mass_,2)) * sin2_2th / (48 * _PI_);
    else
        rate[2] = 0;
        
    if (m >= _charged_pion_mass_ + _electron_mass_)
        rate[3] = pow(_GF_,2) * pow(_pion_decay_,2) * m * sqrt((pow(m,2) - pow(_charged_pion_mass_ + _electron_mass_,2)) * (pow(m,2) - pow(_charged_pion_mass_ - _electron_mass_,2))) * sin2_2th / (16 * _PI_);
    else
        rate[3] = 0;
        
    if (m >= _charged_pion_mass_ + _muon_mass_)
        rate[4] = pow(_GF_,2) * pow(_pion_decay_,2) * m * sqrt((pow(m,2) - pow(_charged_pion_mass_ + _muon_mass_,2)) * (pow(m,2) - pow(_charged_pion_mass_ - _muon_mass_,2))) * sin2_2th / (16 * _PI_);
    else
        rate[4] = 0;
}

// Type 1 is a delta function

double sterile::get_decay_type_one(double energy, double width, double product_mass)
{
    double neutrino_energy = get_monoenergy(m, 0, product_mass);
    double factor = 0;
    double sigma = sigma_delta_boxcar * width;
    double height = 1 / (width - 2 * sigma);
    double sig2 = 1 / (2 * pow(sigma, 2));
    if(energy <= neutrino_energy){
        factor = 0;
    } else if(energy <= neutrino_energy + sigma){
        factor = pow(energy - neutrino_energy, 2) * height * sig2;
    } else if(energy <= neutrino_energy + 2 * sigma){
        factor = height * (1 - pow(energy - neutrino_energy - 2 * sigma, 2) * sig2);
    } else if(energy <= neutrino_energy + width - 2 * sigma){
        factor = height;
    } else if(energy <= neutrino_energy + width - sigma){
        factor = height * (1 - pow(energy - width + 2 * sigma - neutrino_energy, 2) * sig2);
    } else if(energy <= neutrino_energy + width){
        factor = pow(energy - width - neutrino_energy, 2) * height * sig2;
    }

    return factor;    
}

// Type 2 is a boxcar

double sterile::get_decay_type_two(double energy, double spectator_mass)
{
    double gamma_pion;
    double v_pion;
    double p_pion;

    double neutrino_energy_pion = get_monoenergy(_charged_pion_mass_, 0, _muon_mass_);
    double factor = 0;

    compute_kinetics(m, _charged_pion_mass_, spectator_mass, &gamma_pion, &v_pion, &p_pion);
    double width = 2 * gamma_pion * v_pion * neutrino_energy_pion;
    double min = gamma_pion * neutrino_energy_pion * (1 - v_pion);
    double sigma = 0.05 * width;
    double height = 1 / (width - 2 * sigma);
    double sig2 = 1 / (2 * pow(sigma, 2));
    if(energy <= min){
        factor = 0;
    } else if(energy <= min + sigma){
        factor = pow(energy - min, 2) * height * sig2;
    } else if(energy <= min + 2 * sigma){
        factor = height * (1 - pow(energy - min - 2 * sigma, 2) * sig2);
    } else if(energy <= min + width - 2 * sigma){
        factor = height;
    } else if(energy <= min + width - sigma){
        factor = height * (1 - pow(energy - width + 2 * sigma - min, 2) * sig2);
    } else if(energy <= min + width){
        factor = pow(energy - width - min, 2) * height * sig2;
    }

    return factor;
}

double sterile::get_decay_type_three(double energy)
{
    double gamma_muon;
    double v_muon;
    double p_muon;

    double neutrino_max_energy_muon = get_monoenergy(_muon_mass_, 0, _electron_mass_);
    compute_kinetics(m, _muon_mass_, _charged_pion_mass_, &gamma_muon, &v_muon, &p_muon);

    double a = energy / (gamma_muon * (1 + v_muon));
    double b = neutrino_max_energy_muon;

    if (energy / (gamma_muon * (1 - v_muon)) < b)
        b = energy / (gamma_muon * (1 - v_muon));
    
    if (a < b)
        return (get_gamma_a(b) - get_gamma_a(a)) / (2 * gamma_muon * v_muon * (get_gamma_b(neutrino_max_energy_muon) - get_gamma_b(0)));
    else
        return 0.; 
}

double sterile::get_decay_type_four(double energy, double pion_spectator_mass)
{
    double gamma_pion;
    double v_pion;
    double p_pion;
    double gamma_muon;
    double v_muon;
    double p_muon;

    compute_kinetics(m, _charged_pion_mass_, pion_spectator_mass, &gamma_pion, &v_pion, &p_pion);
    compute_kinetics(_charged_pion_mass_, _muon_mass_, 0, &gamma_muon, &v_muon, &p_muon);
    
    double start = gamma_pion * (gamma_muon * _muon_mass_ - v_pion * p_muon);
    double end =  gamma_pion * (gamma_muon * _muon_mass_ + v_pion * p_muon);
    dep_vars* four_vals = new dep_vars(100);
    gel_dummy_vars* x = new gel_dummy_vars(100, start, end);

    for (int i = 0; i < 100; i++){
        four_vals->set_value(i,type_four_integrand(energy, x->get_value(i)));
    }

    double integral = x->integrate(four_vals);
    return integral / (2 * gamma_pion * v_pion * p_muon);
}

void sterile::compute_dPdtdE(dummy_vars* energies_cm, double temp_cm, dep_vars** p_all){
    int num_bins = energies_cm->get_length();

    double* energy_bins = new double[num_bins];
    double* bin_widths = new double[num_bins-1];
    
    for(int i = 0; i < num_bins; i++){
        energy_bins[i] = energies_cm->get_value(i) * temp_cm;
        if (i > 0)
            bin_widths[i-1] = energy_bins[i] - energy_bins[i-1];
                }
    
    double energy, bin_w, d1, d2, d3_2, d3_4, d4_2, d4_3, d4_4;
    
    for(int i = 0; i < num_bins; i++){
        energy = energy_bins[i];
        bin_w = d1 = d2 = d3_2 = d3_4 = d4_2 = d4_3 = d4_4 = 0;
        
        if(i < num_bins-1)
            bin_w = bin_widths[i];
            
        d1 = rate[1] * get_decay_type_one(energy, bin_w, 0);
        
        if (rate[2] > 0)
            d2 = rate[2] * get_decay_type_one(energy, bin_w, _neutral_pion_mass_);
            
        if (rate[3] > 0){
            d3_2 = rate[3] * get_decay_type_two(energy, _electron_mass_);
            d3_4 = rate[3] * get_decay_type_four(energy, _electron_mass_);
        }
        
        if (rate[4] > 0){
            d4_2 = rate[4] * get_decay_type_two(energy, _muon_mass_);
            d4_3 = rate[4] * get_decay_type_three(energy);
            d4_4 = rate[4] * get_decay_type_four(energy, _muon_mass_);
        }
        
        p_all[NU_E]->set_value(i, d1 + d2 + d3_4 + d4_3 + d4_4);
        p_all[NUBAR_E]->set_value(i, d3_4 + d4_3 + d4_4);
        p_all[NU_MU]->set_value(i, d1 + d2 + d3_2 + d3_4 + d4_2 + d4_3 + d4_4);
        p_all[NUBAR_MU]->set_value(i, d3_2 + d3_4 + d4_2 + d4_3 + d4_4);
        p_all[NU_TAU]->set_value(i, d1 + d2);
        p_all[NUBAR_TAU]->set_value(i, 0);
                    
    }
    
    
    delete[] energy_bins;
    delete[] bin_widths;
    return;
}

void sterile::compute_full_term(dummy_vars* energies_cm, double temp_cm, double n_s, double dtda, dep_vars** p_all)
{
    compute_dPdtdE(energies_cm, temp_cm, p_all);
    
    double c = 2 * _PI_ * _PI_ / (temp_cm * temp_cm) * n_s * dtda;
    dep_vars* coeff = new dep_vars(energies_cm->get_length());
    coeff->set_value(0, 0.);
    for(int i = 1; i < energies_cm->get_length(); i++)
        coeff->set_value(i, c / pow(energies_cm->get_value(i), 2));
        
    for(int i = 0; i < 6; i++)
        p_all[i]->multiply_by(coeff);
}


double get_monoenergy(double m0, double m1, double m2){
    return (pow(m0,2) + pow(m1,2) - pow(m2,2)) / (2 * m0);
}

void compute_kinetics(double m0, double m1, double m2, double* gamma, double* v, double* p){
    double particle_energy = get_monoenergy(m0,m1,m2);

    *gamma = particle_energy / m1;
    *p = sqrt(pow(particle_energy,2) - pow(m1,2));
    *v = *(p) / particle_energy;
}

double get_gamma_a(double evaluate_at){
    double term_1 = 3 * pow(_electron_mass_,4);
    double term_2 = 6 * pow(_electron_mass_,2) * _muon_mass_ * evaluate_at;
    double term_3 = pow(_muon_mass_,2) * evaluate_at * (4 * evaluate_at - 3 * _muon_mass_);
    double term_4 = pow(_electron_mass_,4) * _muon_mass_ * log(_muon_mass_ - 2 * evaluate_at);
    return -((evaluate_at / 6) * (term_1 + term_2 + term_3) + term_4 / 4) / pow(_muon_mass_,3);
}

double get_gamma_b(double evaluate_at){
    double term_1 = -16 * pow(_electron_mass_,2) * _muon_mass_ * pow(evaluate_at,3);
    double term_2 = -6 * pow(_electron_mass_,4) * evaluate_at * (_muon_mass_ + evaluate_at);
    double term_3 = -4 * pow(_muon_mass_,2) * pow(evaluate_at,3) * (3 * evaluate_at - 2 * _muon_mass_);
    double term_4 = -3 * pow(_electron_mass_,4) * pow(_muon_mass_,2) * log(_muon_mass_ - 2 * evaluate_at);
    return (term_1 + term_2 + term_3 + term_4) / (24 *pow(_muon_mass_,3));
}

double type_four_integrand(double energy, double muon_energy){
    double gamma_muon = muon_energy / _muon_mass_;
    double v_muon = sqrt(pow(muon_energy,2) - pow(_muon_mass_,2)) / muon_energy;
    double neutrino_max_energy_muon = get_monoenergy(_muon_mass_, 0, _electron_mass_);
    double a = energy / (gamma_muon * (1 + v_muon));
    double b = neutrino_max_energy_muon;

    if (energy / (gamma_muon * (1 - v_muon)) < b){
        b = energy / (gamma_muon * (1 - v_muon));
    }
    if (b - a > 1e-4){
        return (get_gamma_a(b) - get_gamma_a(a)) / (2 * gamma_muon * v_muon * (get_gamma_b(neutrino_max_energy_muon) - get_gamma_b(0)));
    } else {
        return 0;
    }
}





double sterile::energy(double T)
{   return 0;}

double sterile::pressure(double T)
{   return 0;}

double sterile::drho_dT(double T)
{   return 0;}

double sterile::dP_dT(double T)
{   return 0;}