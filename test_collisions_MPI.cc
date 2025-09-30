#include "collisions.hh"
#include "mixed_dummy_vars.hh"
#include "dummy_dep_vars.hh"
#include "freqs.hh"
#include "collisions_MPI.hh"
#include "sterile_decay.hh"

#include "mpi.h"


#include <iostream>
#include <fstream>
#include <chrono>

using std::cout;
using std::endl;

using namespace std;

using namespace std::chrono;

int main(int argc, char* argv[])
{
    //------------- Stuff for MPI -------------
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    //-----------------------------------------
    
    // Stuff for sterile_decay 
    int N = 200;
    double m_s = 300.0;
    
    double a0 = 0.10;
    double af = 0.11;
    double T_0 = 1/a0;
    
    //sterile* nu_s = new sterile(m_s, 5*pow(10,-5));
    
    mixed_dummy_vars* eps = new mixed_dummy_vars(a0, af, m_s, N);
    sterile_decay* sim = new sterile_decay(m_s, 5*pow(10,-5), T_0, eps);
    // ----------------
    
    collisions* C_MPI = new collisions(myid, numprocs, eps);
    
    freqs_ntT* C_values0 = new freqs_ntT(eps, 0., 0., 10., false); // are the inputs for this function okay?
    freqs_ntT* C_valuesf = new freqs_ntT(eps, 0., 0., 10., false); // are the inputs for this function okay?
    
    // Calculate collision integrals before run
    auto start = high_resolution_clock::now();
    
    freqs_ntT* y_vals0 = sim->get_y_values();
    
    C_MPI->C(T_0, y_vals0, true, C_values0); 
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    if(myid == 0){
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
        
        // print C values to a file
        ofstream MyFile1("outputs/C_values0.csv");
        C_values0->print_eps_nus(MyFile1);
        MyFile1.close();
    }
    
    // Run sterile_decay
    sim->run(1500, 1, af, "outputs/test_MPI_output", true);
    
    // Calculate collision integrals after run
    auto start2 = high_resolution_clock::now();
    
    freqs_ntT* y_valsf = sim->get_y_values(); 
    
    C_MPI->C(T_0, y_valsf, true, C_valuesf); 
        
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(stop2 - start2);
    
    if(myid == 0){
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
        
        // print C values to a file
        ofstream MyFile2("outputs/C_valuesf.csv");
        C_valuesf->print_eps_nus(MyFile2);
        MyFile2.close();

    }
    
    delete C_values0;
    delete C_valuesf;
    delete C_MPI;
        
    delete eps;
    
    MPI_Finalize();
    return 0;
}