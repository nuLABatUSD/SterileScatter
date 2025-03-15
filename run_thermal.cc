#include <iostream>
#include <fstream>

#include "sterile_decay.hh"
#include "dummy_dep_vars.hh"
#include "ODESolve.hh"

using std::cout;
using std::endl;
using std::ofstream;

int main()
{
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    double T0 = 20.;
    double a_f = 1./10.;
    
    sterile_decay* standard = new sterile_decay(0., 0., T0, eps);
    standard->run(1000, 1, a_f, "standard.csv", true);
    
    
    sterile_decay* decays = new sterile_decay(300., 1.e-5, T0, eps);
    //decays->run(1000, 1, a_f, "decays.csv", true);
    
    gel_linspace_gl* eps_prime = new gel_linspace_gl(2., 20., 200);
    
    sterile_decay* decay_prime = new sterile_decay(300., 1.e-5, T0, eps_prime);
//    decay_prime->run(1000, 1, a_f, "decay_prime.csv", true);

    sterile_decay* dec2 = new sterile_decay(300., 5.e-5, 1./T0, a_f, 200);
    
    
    dec2->run(2000, 2, a_f, "decay-0.csv", true);
    
    dec2->shift_eps(1./5.);
    dec2->run(2000, 2, 0.2, "decay-1.csv", true);
    
    ofstream file;
    file.open("dec_eps.csv");
    dec2->print_eps_file(file);
    file.close();    
    
    delete dec2;
    
    delete eps_prime;
    delete decay_prime;
    
    delete standard;
    delete decays;
    delete eps;
    return 0;
}