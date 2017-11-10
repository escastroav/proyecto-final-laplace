
#include <iostream>
#include <vector>
#include <cmath>


// constants
const double DELTA = 0.1;
const double L = 2.0; 
const int N = int(L/DELTA)+1;
const int STEPS = 200;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void evolve(Matrix & m);
void print(const Matrix & m);
void init_gnuplot(void);
void plot_gnuplot(const Matrix & m);


