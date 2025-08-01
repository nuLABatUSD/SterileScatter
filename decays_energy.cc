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

int main()
{
    string folder_name = "decaysANDenergy";
    int result = mkdir(folder_name.c_str(), 0777);
    
    double m_s = 300.0;
    double theta = 5e-5;
    int N = 200;
    
    double a0 = 1.0; // is this 1, since Tcm will be 1?
    double af = 1.01; // what do we set af to?
    double Tcm = 1/a0;
    
    ofstream MyFile1(folder_name + "/eps_output.csv"); // Output file for eps & weights of each run
    ofstream MyFile2(folder_name + "/decays_output.csv"); // Output file for neutrino energy densities pre & post shift + ratio
    
    
    sterile* nu_s = new sterile(m_s, theta);
    mixed_dummy_vars* eps = new mixed_dummy_vars(a0,af,m_s,N);
    
    // create pointer to pointer dep_vars**
    dep_vars** p_all = new dep_vars*[6];
    for (int i = 0; i < 6; i++)
        p_all[i] = new dep_vars(N);
    
    // compute dPdtdE
    nu_s->compute_dPdtdE(eps, Tcm, p_all);
    
    // print p_all to file
    for (int j = 0; j < N; j++){
        for (int i = 0; i < 6-1; i++){
            MyFile2 << p_all[i]->get_value(j) << ",";
        }
        MyFile2 << p_all[6-1]->get_value(j) << endl;
    }
    
    // print eps to file
    eps->print_csv(MyFile1);
    
    // delete dep_vars**
    for(int i = 0; i < 6; i++)
        delete p_all[i];
    delete[] p_all;
    
    MyFile1.close();
    MyFile2.close();
    
    return 1;    
}
