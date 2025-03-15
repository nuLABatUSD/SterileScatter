#ifndef __DUMMY_DEP_HH__
#define __DUMMY_DEP_HH__

#include <iostream>

using std::ostream;

class dep_vars;

class dummy_vars{
    protected:
        int N;
        double* values;
        double* weights;
        
        double max_linspace;
    public:
        dummy_vars();
        dummy_vars(int);
        dummy_vars(dummy_vars*);
        ~dummy_vars();
        
        double get_value(int);
        double get_weight(int);
        int get_length();
        double get_max_linspace();
        
        void set_value(int, double);
        void set_weight(int, double);
        
        void set_trap_weights();
        
        int bin_below(double);
        
        double integrate(dep_vars*);
        double partial_integrate_end(int, dep_vars*);
        double integrate_pow(dep_vars*, double);
        double partial_integrate_pow_end(int, dep_vars*, double);
        
        void print_all();
        void print_csv(ostream&);
};

class gl_dummy_vars : public dummy_vars{
    public:
        gl_dummy_vars(int);
        gl_dummy_vars(int, double);
};

class gel_dummy_vars : public dummy_vars{
    public:
        gel_dummy_vars(int, double, double);
};

class linspace_and_gl : public dummy_vars{
    protected:
        int num_lin;
        int num_gl;
        
    public:
        linspace_and_gl(double, double, int, int);
        linspace_and_gl(linspace_and_gl*);
            
        int get_num_lin();
        int get_num_gl();
};

class linspace_for_trap : public linspace_and_gl{
    public:
        linspace_for_trap(double, double, int);
        linspace_for_trap(linspace_for_trap*);
};

class gel_linspace_gl : public dummy_vars{
    protected:
        const int default_N_gel = 10;
        const int default_N_gl = 5;
        const double max_lin_sm = 20.;
        
        int num_gel, num_lin, num_gl;
        
    public:
        gel_linspace_gl(double, double, int);
        gel_linspace_gl(gel_linspace_gl*);
        
        int get_num_gel();
        int get_num_lin();
        int get_num_gl();
        double get_min_linspace();
        double get_delta_linspace();
    
};

class gel_inner_integral : public dummy_vars{
    protected:
        int N_gel;
        int N_bins_partition;
        
    public:
        gel_inner_integral(dummy_vars* eps, double xmax, int Ngel = 5, int Nb = 10);
};


class dep_vars{
    protected:
        int N;
        double* values;
        
    public:
        dep_vars(int);
        dep_vars(double*, int);
        dep_vars(dep_vars*);
        ~dep_vars();
        
        double get_value(int);
        int get_length();
        
        bool isnan();
        
        void set_value(int, double);
        
        void zeros();
        void copy(dep_vars*);

        void multiply_by(double);
        void multiply_by(dep_vars*);
        void add_to(double, dep_vars*);
        
        void print_all();
        void print(int N_top = 3, int N_bot = 1);
        void print_csv(ostream& os);
};
#endif