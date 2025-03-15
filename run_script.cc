#include <iostream>
#include <fstream>
#include <string>
#include <chrono>


#include "sterile_decay.hh"
#include "dummy_dep_vars.hh"
#include "ODESolve.hh"

using std::cout;
using std::endl;
using std::ofstream;
using namespace std::chrono;

int main(int argc, char* argv[])
{
    double T0 = 10.;
    double a_mult = 1.5;
    double a_end = 100.;

    double af = a_mult / T0;
    
    std::string folder_name(argv[1]);
    
    std::string eps_file_name = folder_name + "/decay-eps.csv";
    std::string shift_file_name = folder_name + "/decay-shift.csv";
        
    ofstream eps_file;
    eps_file.open(eps_file_name);
    
    ofstream shift_file;
    shift_file.open(shift_file_name);
    
    double mass = std::atof(argv[2]);
    double th = std::atof(argv[3]);
    
    sterile_decay* dec = new sterile_decay(mass, th, 1./T0, af, 200);
    
    dec->print_csv(shift_file);
    dec->print_eps_file(eps_file);
    
    auto start = high_resolution_clock::now();

    for(int i = 0; i < 50; i++){
        std::string file_name = folder_name + "/decay-" + std::to_string(i) + ".csv";
        dec->run(2000, 10, af, file_name, true);        
        dec->print_csv(shift_file);

        af = dec->shift_eps_by_multiple(a_mult);
        
        dec->print_eps_file(eps_file);
        dec->print_csv(shift_file);
        
        if (af > a_end * a_mult)
            break;
    }
    
    cout << endl << "N_eff = " << dec->get_Neff() << endl << endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << endl << "Time elapsed: "
     << duration.count()/1000./60. << " minutes" << endl;
    
    eps_file.close();
    shift_file.close();
    
    delete dec;
    return 0;
}