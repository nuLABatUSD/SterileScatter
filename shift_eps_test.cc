#include <iostream>
#include <cmath>
#include <fstream>
#include "dummy_dep_vars.hh"
#include "mixed_dummy_vars.hh"
#include "sterile.hh"
#include "sterile_decay.hh"
#include "ODESolve.hh"

#include <string>
using std::string;
using std::to_string;

int main()
{
    /*
    // -------------- Simulation ------------------
    int N = 200;
    double m_s = 300.0;
    
    sterile* nu_s = new sterile(m_s, 5*pow(10,-5));
    
    double a0 = 0.1;
    double af = 0.11;
    mixed_dummy_vars* eps = new mixed_dummy_vars(0,a0,af,m_s,N);
    
    
    double T_0 = 1/a0;
    double T_f = 1/af;
    sterile_decay* sim = new sterile_decay(m_s, 5*pow(10,-5), T_0, eps);
    
    string sim_output = "outputs/sim_output" + to_string(1) + ".csv";
    string eps_output = "outputs/eps_output" + to_string(1) + ".csv";
    
    sim->run(1500, 1, af, sim_output, true);
    ofstream MyFile1(eps_output);
    sim->print_eps_file(MyFile1);
    
    MyFile1.close();

    string sim_output2 = "outputs/sim_output" + to_string(2) + ".csv";
    string eps_output2 = "outputs/eps_output" + to_string(2) + ".csv";
    
    double new_af = sim->shift_eps(af + 0.01);
    
    sim->run(1500, 1, new_af, sim_output2, true);
    ofstream MyFile2(eps_output2);
    sim->print_eps_file(MyFile2);
    
    MyFile2.close();
    */
    
    int N = 200;
    double m_s = 300.0;
    
    sterile* nu_s = new sterile(m_s, 5*pow(10,-5));
    
    double a0, af, T_0, T_f, e1, e2, ratio;
    
    a0 = 0.1;
    af = 0.11;
    T_0 = 1/a0;
    T_f = 1/af;

    mixed_dummy_vars* eps = new mixed_dummy_vars(0,a0,af,m_s,N);
        
    sterile_decay* sim = new sterile_decay(m_s, 5*pow(10,-5), T_0, eps);

    ofstream MyFile("outputs/eps_output.csv"); // Output file for eps & weights of each run
    ofstream MyFile2("outputs/ned_output.csv"); // Output file for neutrino energy densities pre & post shift + ratio
    
    int num_runs = 100;
    for(int i = 0; i < num_runs; i++){
        a0 = 0.1 + (0.01 * i);
        af = 0.11 + (0.01 * i);
        T_0 = 1/a0;
        T_f = 1/af;
    
        string sim_output = "outputs/sim_output" + to_string(i+1) + ".csv";
        
        sim->run(1500, 1, af, sim_output, true);
        
        
        sim->print_eps_file(MyFile);
        
        
        e1 = sim->get_neutrino_energy(af);
        MyFile2 << e1 << ",";
        
        sim->shift_eps_by_multiple(1.1);
        
        e2 = sim->get_neutrino_energy(af);
        MyFile2 << e2 << ",";
        
        ratio = e1/e2;
        MyFile2 << ratio << endl;
    }
    
    MyFile.close();
    MyFile2.close();
    
    
    
    return 1;    
}
