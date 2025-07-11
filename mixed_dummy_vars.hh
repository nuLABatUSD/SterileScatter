#ifndef __MIXED_DUMMY_VARS_HH__
#define __MIXED_DUMMY_VARS_HH__

#include "dummy_dep_vars.hh"

#include <iostream>

using std::ostream;

// ----------- Albert Edits: 6/24/25 ---------------
class mixed_dummy_vars : public dummy_vars{
    protected:
    const int default_N_gl = 5; // 5 Gauss-Laguerre points for the highest eps values
    
    const int def_N_ints = 4; // number of G-Leg + LS intervals

    double Estart;
    double a0;
    double af;
    double ms;
    int N;
    
    double key_energies[6];
    
    public:
    mixed_dummy_vars(double, double, double, double, int);// E_start, a_0, a_f, m_s, number of points num
    mixed_dummy_vars(mixed_dummy_vars*);
    
    double get_ms();
    double get_key_energies(int);// method to get the 6 important energies, need to input index of energy you want?

};
// ----------- End of Albert Edits -----------------

#endif