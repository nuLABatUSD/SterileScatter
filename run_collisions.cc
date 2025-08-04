#include "collisions.hh"
#include "mixed_dummy_vars.hh"
#include "dummy_dep_vars.hh"
#include "freqs.hh"

#include <iostream>
#include <chrono>

using std::cout;
using std::endl;

using namespace std::chrono;

int main()
{
    mixed_dummy_vars* eps = new mixed_dummy_vars(0.10, 0.11, 300., 200);
    nu_nu_collision* C;
    
    double results[6];
    
    auto start = high_resolution_clock::now();
    
    for(int i = 0; i < eps->get_length(); i++){
        C = new nu_nu_collision(i, eps);
        C->compute_R(10., 10., results);
        
        cout << i << ", eps = " << eps->get_value(i) << ", R = " << results[0] << endl;        
        
        delete C;
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
        
    delete eps;
}