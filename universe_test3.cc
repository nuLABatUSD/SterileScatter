#include <iostream>
#include <fstream>
#include "dummy_dep_vars.hh"
#include "sterile.hh"
#include "sterile_decay.hh"
#include "ODESolve.hh"


//    Following Chad's Notes p.4-5

int main()
{
    linspace_and_gl* eps = new linspace_and_gl(0,20,201,5);
        // dummy vars w/ 201 linspace 0->20 5 Gauss-Laguerre pts
        
    double T_0 = 20.0;
    sterile_decay* standard = new sterile_decay(0, 0, T_0, eps);
    
    standard->run(1000, 1, 100.0, "output3.csv", true);
        // up to 1000 print file lines, 1 each step, final:a_f, output filename, verbose... print things to terminal
    
    ofstream MyFile("eps_output.csv");
    standard->print_eps_file(MyFile); // prints eps to file
   
   MyFile.close();
   delete eps;
   delete standard;
    

    return 1;    
}

