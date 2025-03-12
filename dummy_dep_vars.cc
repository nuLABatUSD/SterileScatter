#include <iostream>
#include <cmath>
#include <iomanip>

#include "dummy_dep_vars.hh"
#include "gl_vals.hh"
#include "gel_vals.hh"

using std::cout;
using std::endl;
using std::ostream;

//dummy_vars
dummy_vars::dummy_vars()
{   N = -1; }

dummy_vars::dummy_vars(int num){
    N = num;
    values = new double[N]();
    weights = new double[N]();
}

dummy_vars::dummy_vars(dummy_vars* copy_me)
{
    N = copy_me->get_length();
    values = new double[N]();
    weights = new double[N]();

    for(int i = 0; i<N; i++)
        {
            values[i] = copy_me->get_value(i);
            weights[i] = copy_me->get_weight(i);
        }
}

dummy_vars::~dummy_vars(){
    delete[] values;
    delete[] weights;
}

double dummy_vars::get_value(int i){
    return values[i];
}

double dummy_vars::get_weight(int i)
{ return weights[i]; }

int dummy_vars::get_length(){
    return N;
}

void dummy_vars::set_value(int i, double v)
{values[i] = v;}

void dummy_vars::set_weight(int i, double w)
{weights[i] = w;}

void dummy_vars::set_trap_weights(){
   weights[0] = 0.5 * (values[1] - values[0]);
   weights[N-1] = 0.5 * (values[N-1] - values[N-2]);
    for(int i=1; i<N-1; i++){
       weights[i] = 0.5 * (values[i+1] - values[i-1]);
   }
}

int dummy_vars::bin_below(double x){
    if(N < 2)
        return N;
    int guess = (int) ( (x - values[0]) / ((values[N-1] - values[0]) / (N-1)) );
        
    if (guess >= N-1){
        if (values[N-1] <= x)
            return N-1;
        guess = N-2;
    }
    else if( values[guess] <= x && x < values[guess+1])
        return guess;
        
    int pm = 1;
    if(values[guess] > x)
        pm = -1;
        
    int g;
    for(int i = 0; i < N; i++){
        g = guess + pm * i;
        if (g == 0 || g == N-1)
            return g;
            
        if (values[g] <= x && x < values[g+1])
            return g;
    }
    return g;
    
}

double dummy_vars::integrate(dep_vars* fvals){
    double result = 0;
    for (int i = 0; i<N; i++){
       result += fvals->get_value(i) * weights[i]; 
    }
    return result;    
}

void dummy_vars::print_all(){
    for(int i =0; i<N; i++){
        cout << values[i] << endl;
    }
}

gl_dummy_vars::gl_dummy_vars(int num_gl) : gl_dummy_vars(num_gl, 0.)
{ ;
/*    switch(num_gl){
        case 2:
            for(int i=0; i<N; i++){
                values[i] = xvals_2[i];
                weights[i] = wvals_2[i] * exp(xvals_2[i]);
             }
            break;
        case 5:
            for(int i=0; i<N; i++){
                values[i] = xvals_5[i];
                weights[i] = wvals_5[i] * exp(xvals_5[i]);
             }
            break;
        case 10:
             for(int i=0; i<N; i++){
                values[i] = xvals_10[i];
                weights[i] = wvals_10[i] * exp(xvals_10[i]);
             }
            break;
        case 50:
             for(int i=0; i<N; i++){
                values[i] = xvals_50[i];
                weights[i] = wvals_50[i] * exp(xvals_50[i]);
             }
            break;
        default:
            cout << "Error: this Gauss Laguerre number is not supported" << endl;
                
    }
    */
}

gl_dummy_vars::gl_dummy_vars(int num_gl, double start) : dummy_vars(num_gl)
{
    const double* val;
    const double* w;
        
    switch(num_gl){
        case 2:
            val = xvals_2;
            w = wvals_2;
            break;
        case 5:
            val = xvals_5;
            w = wvals_5;
            break;
        case 10:
            val = xvals_10;
            w = wvals_10;
            break;
        case 50:
            val = xvals_50;
            w = wvals_50;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;    
    }
    
    for(int i = 0; i < num_gl; i++){
        values[i] = val[i] + start;
        weights[i] = w[i] * exp(val[i]);
    }
}

gel_dummy_vars::gel_dummy_vars(int num_gel, double start, double end) : dummy_vars(num_gel)
{
    double half_width = (end - start) / 2.;
    double slope_shift = half_width;
    double shift = (end + start) / 2.;
    
    const double* val;
    const double* w;
        
    switch(num_gel){
        case 2:
            val = gel_vals_2;
            w = gel_weights_2;
            break;
        case 5:
            val = gel_vals_5;
            w = gel_weights_5;
            break;
        case 10:
            val = gel_vals_10;
            w = gel_weights_10;
            break;
        case 50:
            val = gel_vals_50;
            w = gel_weights_50;
            break;
        case 100:
            val = gel_vals;
            w = gel_weights;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;
    }
    
    for(int i = 0; i < num_gel; i++){
        values[i] = slope_shift * val[i] + shift;
        weights[i] = half_width * w[i];
    }
    
}

linspace_and_gl::linspace_and_gl(double xmin, double xmax, int numlin, int numgl) : dummy_vars(numlin+numgl)
{
    num_lin = numlin;
    num_gl = numgl;
    
    double dx_val = (xmax - xmin) / (num_lin -1);
    for (int i = 0; i<num_lin; i++){
        values[i] = xmin + dx_val * i;
        weights[i] = dx_val;
    }
    
    weights[0] = dx_val / 2;
    weights[num_lin-1] = dx_val / 2;

    if (numgl > 0){
        gl_dummy_vars* gl = new gl_dummy_vars(numgl, xmax);
        
        for (int i = 0; i < numgl; i++){
            values[num_lin+i] = gl->get_value(i);
            weights[num_lin+i] = gl->get_weight(i);
        }
    }
}

int linspace_and_gl::get_num_lin()
{   return num_lin; }

int linspace_and_gl::get_num_gl()
{   return num_gl;  }

double linspace_and_gl::get_max_linspace()
{   return values[num_lin-1]; }

linspace_and_gl::linspace_and_gl(linspace_and_gl* copy_me) : dummy_vars(copy_me)
{
    num_lin = copy_me->get_num_lin();
    num_gl = copy_me->get_num_gl();
}

linspace_for_trap::linspace_for_trap(double xmin, double xmax, int num) : linspace_and_gl(xmin, xmax, num, 0)
{ ; }

linspace_for_trap::linspace_for_trap(linspace_for_trap* copy_me) : linspace_and_gl(copy_me)
{ ;}

gel_inner_integral::gel_inner_integral(dummy_vars* eps, double xmax, int Ngel, int Nb) : dummy_vars() {
    gel_dummy_vars* gel = new gel_dummy_vars(Ngel, 0., 1.);
    if (gel->get_weight(0) == 0){
        cout << "gel_inner_integral failed to create Gauss-Legendre dummy_vars with N = " << N << endl;
        delete gel;
        return;
    }
    delete gel;
    
    int ind = eps->bin_below(xmax);
    
    int Ndv = ind / Nb;
    if (ind % Nb >= Nb / 2)
        Ndv++;
        
    N = Ndv * Ngel;
    values = new double[N]();
    weights = new double[N]();
    
    for(int i = 0; i < Ndv; i++){
        if(i < Ndv - 1)
            gel = new gel_dummy_vars(Ngel, eps->get_value(i*Ngel), eps->get_value((i+1)*Ngel));
        else
            gel = new gel_dummy_vars(Ngel, eps->get_value(i*Ngel), xmax);
        for(int j = 0; j < Ngel; j++){
            values[i*Ngel + j] = gel->get_value(j);
            weights[i*Ngel + j] = gel->get_weight(j);
        }
        delete gel;
    }
}

dep_vars::dep_vars(int size)
{
    N = size;
    values = new double[N]();
}
    
dep_vars::dep_vars(double* copy_me, int size)
{
    N = size;
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me[i];
}
    
dep_vars::dep_vars(dep_vars* copy_me)
{
    N = copy_me->get_length();
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me->get_value(i);
}
    
dep_vars::~dep_vars()
{delete[] values;}

double dep_vars::get_value(int i)
{return values[i];}

int dep_vars::get_length()
{return N;}

void dep_vars::set_value(int i, double v)
{values[i] = v;}

void dep_vars::zeros()
{
    for (int i = 0; i < N; i++)
        values[i] = 0.0;
}

void dep_vars::copy(dep_vars* z)
{
    for (int i = 0; i < N; i++) 
        values[i] = z -> get_value(i);
}

void dep_vars::multiply_by(double scalar)
{
    for (int i = 0; i < N; ++i) 
        values[i] *= scalar;
}

void dep_vars::multiply_by(dep_vars* vector)
{
    if(vector->get_length() != N){
        cout << "Error using dep_vars::multiply_by(dep_vars*). Two dep_vars objects have different lengths." << endl;
        return;
    }
    
    for(int i = 0; i < N; i++)
        values[i] *= vector->get_value(i);
}

void dep_vars::add_to(double c, dep_vars* z)
{
    for (int i = 0; i < N; ++i)
        values[i] += c * z -> get_value(i);
}

bool dep_vars::isnan(){
    for(int i = 0; i < N; i++)
        if(std::isnan(values[i])){
            cout << "ERROR: dep_vars object is nan at index " << i << endl;
            return true;
        }
    return false;
}

void dep_vars::print_all()
{
    for (int i = 0; i < N; i++)
        cout << values[i] << endl;
}
    
void dep_vars::print(int N_top, int N_bot)
{
    if (N <= N_top + N_bot)
        print_all();
    else if (N_top < 0 || N_bot < 0)
        print_all();
    else
    {
        for (int i = 0; i < N_top; i++)
            cout << values[i] << endl;
        cout << "..." << endl;
        for (int i = 0; i < N_bot; i++)
            cout << values[N - N_bot + i] << endl;
    }
}

void dep_vars::print_csv(ostream& os)
{
    for (int i = 0; i < N-1; i++)
        os << values[i] << ", ";
    os << values[N-1];
}
