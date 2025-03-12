#include <iostream>

#include "sterile_decay.hh"
#include "dummy_dep_vars.hh"
#include "ODESolve.hh"

using std::cout;
using std::endl;

int main()
{
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    double T0 = 20.;
    double a_f = 1./10.;
    
    sterile_decay* standard = new sterile_decay(0., 0., T0, eps);
    standard->run(1000, 1, a_f, "standard.csv", true);
    
    
    sterile_decay* decays = new sterile_decay(300., 1.e-5, T0, eps);
    decays->run(1000, 1, a_f, "decays.csv", true);
    
    delete standard;
    delete decays;
    delete eps;
    return 0;
}