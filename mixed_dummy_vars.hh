#ifndef __MIXED_DUMMY_VARS_HH__
#define __MIXED_DUMMY_VARS_HH__

#include "dummy_dep_vars.hh"

#include <iostream>

using std::ostream;

// ----------- Albert Edits: 6/24/25 ---------------
class mixed_dummy_vars : public dummy_vars{
    protected:
    double a0;
    double af;
    double ms;
    int N;
    
    int def_N_ints = 4; // number of G-Leg + LS intervals
    int def_N_gel1 = 10; // default # of points for first gel interval
    int def_N_gel = 5; // default # of points for later gel intervals
    int def_N_linTop = 50; // default # of points for linspace covering Top Hats
    int def_N_linDel = 35; // default # of points for linspace covering Delta Functions
    
    
    double key_energies[6];
    double intBounds[11];
    //double intBounds_SHIFT[11]; // array with 11 points to have additional linspace regions after a shift
    
    public:
    mixed_dummy_vars(double, double, double, int);//a_0, a_f, m_s, number of points num
    
    // This constructor is meant to be called when shifting the eps
    mixed_dummy_vars(double, double, double, double, int); //a_0, a_f, m_s, num, shift = T/F
    
    mixed_dummy_vars(mixed_dummy_vars*);
    
    void calc_energies();
    
    double get_ms();
    double get_a0();
    double get_key_energies(int);// method to get the 6 important energies, need to input index of energy you want?
    double get_intBounds(int); // method to get intBounds
    
    

};
// ----------- End of Albert Edits -----------------




#endif