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

    
    public:
    mixed_dummy_vars(double, double, double, double, int); 
    // E_start, Tcm_0, Tcm_f, m_s, number of points N
    

};
// ----------- End of Albert Edits -----------------

#endif