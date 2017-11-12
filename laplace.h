#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "papi.h"




// constants
const double DELTA = 0.01;
const double L = 10.0; 
const int N = int(L/DELTA)+1;
const int STEPS = 100;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void evolve(Matrix & m);
void print(const Matrix & m);
void init_gnuplot(void);
void plot_gnuplot(const Matrix & m);
void PAPI_measure_Init(int retval, float & ireal_time, float & iproc_time, float & imflops, long long & iflpops);
void PAPI_measure_Final(int retval, float & real_time, float & proc_time, float & mflops, long long & flpops);


