#include <iostream>

#include "dummy_dep_vars.hh"
#include "universe.hh"
#include "sterile.hh"

using std::cout;
using std::endl;

int main()
{
    sterile* ms = new sterile(300., 2.e-5);
    cout << "Sterile mass " << ms->mass() << ", lifetime = " << ms->get_lifetime_s() << endl;
    cout << "Rate = " << ms->get_rate() << endl;
    
    linspace_for_trap* eps = new linspace_for_trap(0., 151., 1000);
    
    dep_vars** dPdtdE = new dep_vars*[6];
    for (int i = 0; i < 6; i++)
        dPdtdE[i] = new dep_vars(eps->get_length());
        
    ms->compute_dPdtdE(eps, 1., dPdtdE);
    
    for (int i = 0; i < 6; i++)
        cout << eps->integrate(dPdtdE[i]) << endl;
        
    delete ms;
    delete eps;
    for (int i = 0; i < 6; i++)
        delete dPdtdE[i];
    delete[] dPdtdE;
    return 0;
}