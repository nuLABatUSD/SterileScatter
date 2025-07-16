#include <cmath>

#include "freqs.hh"
#include "dummy_dep_vars.hh"
#include "constants.hh"

freqs_ntT::freqs_ntT(freqs_ntT* copy_me) : dep_vars(copy_me){
    num_bins = copy_me->get_num_bins();
    eps = new mixed_dummy_vars(copy_me->get_eps()); 
}

freqs_ntT::freqs_ntT(mixed_dummy_vars* e, double n0, double time, double Temp, bool thermal_dist) : dep_vars(6 * e->get_length() + 3){
    num_bins = e->get_length();
    eps = new mixed_dummy_vars(e);
    
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

mixed_dummy_vars* freqs_ntT::get_eps()
{   return eps; }

void freqs_ntT::new_eps(mixed_dummy_vars* e)
{
    delete eps;
    eps = new mixed_dummy_vars(e);
}

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

double freqs_ntT::get_eps_value(int i)
{   return eps->get_value(i);  }

void freqs_ntT::set_f_value(int bin, int type, double val)
{   values[type * num_bins + bin] = val; }

void freqs_ntT::set_sterile_density(double n0)
{   values[6 * num_bins] = n0; }

void freqs_ntT::set_time(double t0)
{   values[6 * num_bins + 1] = t0; }

void freqs_ntT::set_Temp(double T0)
{   values[6 * num_bins + 2] = T0; }

void freqs_ntT::print_eps_file(ostream& os)
{   eps->print_csv(os);  
    os << std::endl; }

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

void freqs_ntT::get_neutrino_distribution(int type, dep_vars* f_out){
    for(int i = 0; i < num_bins; i++)
        f_out->set_value(i, get_f_value(i, type));
}

void freqs_ntT::interpolated_f_values(double x, double* results){
    interpolated_f_values(x, eps->bin_below(x), results);
}

void freqs_ntT::interpolated_f_values(double x, int bin, double* results)
{
    if(x < eps->get_max_linspace()){
        /*int ind = std::max(0, bin-2);
        ind = std::min(num_bins-1-4, ind);*/
        int ind = bin;
        
        double x_vals[5];
        double y_vals[5];
        
        for(int i = 0; i < 5; i++)
            x_vals[i] = eps->get_value(ind+i);
        
        for(int j = 0; j < 6; j++){
            for(int i = 0; i < 5; i++)
                y_vals[i] = get_f_value(ind+i, j);
            results[j] = interpolate_log_fifth(x, x_vals, y_vals);
        }    
    }
    else{
        for(int j = 0; j < 6; j++){
            if(bin == num_bins-1)
                results[j] = interpolate_log_linear(x, eps->get_value(num_bins-2), eps->get_value(num_bins-1), get_f_value(num_bins-2, j), get_f_value(num_bins-1, j));
            else
                results[j] = interpolate_log_linear(x, eps->get_value(bin), eps->get_value(bin+1), get_f_value(bin, j), get_f_value(bin+1, j));
        }
    }
}

void freqs_ntT::interpolate_extrapolate(double x, double T_cm, double* results){
    int bin = eps->bin_below(x);    
    
    int bin0;
    
    double key_eps[6];
    int key_bins[6];
    
    // Calculate key_eps and bin_eps
    for(int i = 0; i < 6; i++){
        key_eps[i] = eps->get_key_energies(i)/T_cm;
        key_bins[i] = eps->bin_below(key_eps[i]);
    }// Improvement: Calculating this every time is silly, these should be stored somewhere and just referenced
    
   /* if(bin == 0){
        for(int i = 0; i < 6; i++){
            std::cout << key_eps[i] << "         ";
        }
        std::cout << std::endl;
        for(int i = 0; i < 6; i++){
            std::cout << key_bins[i] << "         ";
        }
        std::cout << eps->get_max_linspace() << std::endl;
    }*/
    
    //std::cout << "bin_below= " << bin;
    // Check for case: need to check for all 6 key bins
    for(int i = 0; i < 6; i++){
        // Check for case 1: Above key bin
        if(bin == key_bins[i] || bin == key_bins[i]+1 || bin == key_bins[i]+2){
            bin0 = key_bins[i] + 1;
            //std::cout << std::endl << "Case 1 at bin " << bin << std::endl;
            break;
        }
        // Check for case 2: Below or at key bin
        else if(bin == key_bins[i]-1){
            bin0 = key_bins[i] - 4;
            //std::cout << std::endl << "Case 2 at bin " << bin << std::endl;
            break;
        }
        // Neither case is true
        else{
            bin0 = std::max(0,bin-2); // deals w/ x too small
            bin0 = std::min(num_bins-1-4, bin0); // deals w/ x too big
        }
    }
    //std::cout << "     bin0 = " << bin0 << std::endl;
    
    // Call interpolate function
    interpolated_f_values(x, bin0, results);
}

double fifth_order_fit(double x, double* x_vals, double* y_vals){
    double fit = 0;
    double Lj;

    for(int j=0; j<5; j++){
        Lj = 1.0;
        for(int i=0; i<5; i++){
            if(i != j){
                Lj *= (x - x_vals[i]) / (x_vals[j] - x_vals[i]);
            }
        }
        fit += y_vals[j] * Lj;
    }
    return fit;
}

double interpolate_log_fifth(double x, double* x_vals, double* y_vals){
    double y_temp;

    double* y_log = new double[5]();
    for(int j=0; j<5; j++)
        y_log[j] = log(y_vals[j]);
    
    y_temp = fifth_order_fit(x, x_vals, y_log);
    
    delete[] y_log;
    return exp(y_temp);
}

double linear(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0**" << x << std::endl;}
    double slope = (y2-y1)/(x2-x1);
    return slope * (x-x2) + y2;
}

double interpolate_log_linear(double x, double x_val1, double x_val2, double y_val1, double y_val2){
    double y_temp = linear(x, x_val1, x_val2, log(y_val1), log(y_val2));
    return exp(y_temp);
}
