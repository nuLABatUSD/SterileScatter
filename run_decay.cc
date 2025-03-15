#include <iostream>
#include <fstream>
#include <string>

#include "sterile_decay.hh"
#include "dummy_dep_vars.hh"
#include "ODESolve.hh"

using std::cout;
using std::endl;
using std::ofstream;

int main()
{
    double T0 = 20.;
    
    double a_mult = 1.5;
    double af = a_mult / T0;
    
    double a_end = 100.;
    
    ofstream eps_file;
    eps_file.open("decay-eps.csv");
    
    ofstream shift_file;
    shift_file.open("decay-shift.csv");
    
    sterile_decay* dec = new sterile_decay(300., 5.e-5, 1./T0, af, 200);
    
    dec->print_csv(shift_file);
    dec->print_eps_file(eps_file);
    
    for(int i = 0; i < 10; i++){
        std::string file_name = "decay-" + std::to_string(i) + ".csv";
        dec->run(2000, 10, af, file_name, true);        
        dec->print_csv(shift_file);

        dec->shift_eps_by_multiple(a_mult);
        
        dec->print_eps_file(eps_file);
        dec->print_csv(shift_file);

        af *= a_mult;
    }
    
    eps_file.close();
    shift_file.close();
    
    delete dec;
    return 0;
}