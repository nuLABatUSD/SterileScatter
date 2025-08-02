#include <iostream>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"
#include "sterile.hh"
#include "sterile_decay.hh"
#include "ODESolve.hh"

#include <sys/stat.h> // Required for mkdir
#include <sys/types.h> // Required for mode_t

#include <string>
using std::string;
using std::to_string;

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
    // program name, mass, mixing angle, folder name
    
    cout << "Program Name Is: " << argv[0] << endl;
    
    if (argc >= 2){
        cout << "Number of Arguments Passed: " << argc << endl;
        cout << "----Following Are the Command Line Arguments Passed----" << endl;
        cout << "argv[0]: " << argv[0] << endl;
        cout << "Mass - \targv[1]: " << argv[1] << endl;
        cout << "Mixing Angle - \targv[2]: " << argv[2] << endl;
        cout << "Folder Name - \targv[3]: " << argv[3] << endl;
    }
    
    int N = 200;
    double m_s = std::stod(argv[1]);
    double theta = std::stod(argv[2]);
    string folder_name = argv[3];
    
    int result = mkdir(folder_name.c_str(), 0777);
    
    
    sterile* nu_s = new sterile(m_s, theta);
    
    double a0, af, T_0, T_f, e1, e2, ratio;
    
    a0 = 0.1;
    af = a0 + 0.01;
    T_0 = 1/a0;
    T_f = 1/af;

    mixed_dummy_vars* eps = new mixed_dummy_vars(a0,af,m_s,N);
    sterile_decay* sim = new sterile_decay(m_s, theta, T_0, eps);
    

    ofstream MyFile1(folder_name + "/eps_output.csv"); // Output file for eps & weights of each run
    ofstream MyFile2(folder_name + "/Neff_output.csv"); // Output file for neutrino energy densities pre & post shift + ratio
    
    int num_runs = 2000;
    bool decay_state;
    
    for(int i = 0; i < num_runs; i++){
        
    
        string sim_output = folder_name + "/sim_output" + to_string(i+1) + ".csv";
        
        sim->print_eps_file(MyFile1);
        sim->run(1500, 1, af, sim_output, true); // usually run with 1500
        
        
        
        
        sim->shift_eps(af+0.01);
        
        decay_state = sim->check_decay_on();
        
        if(decay_state == false){
            af = 1000;
            
            string sim_output = folder_name + "/sim_output" + to_string(i+1) + ".csv";
            sim->run(5000, 100, af, sim_output, true);
            // need to make sure that this final run is going all the way to af = 1000
            // with 1500, 1, only went to ~24
            
            break;
        }
        
        a0 += 0.01;
        af += 0.01;
    }
    
    MyFile2 << sim->get_Neff() << "," << nu_s->get_lifetime() << "," << nu_s->get_lifetime_s();
    
    MyFile1.close();
    MyFile2.close();
    
    return 1;    
}
