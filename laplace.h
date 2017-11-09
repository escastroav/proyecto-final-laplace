#include <iostream>
#include <vector>
#include <cmath>

// constants
const double DELTA = 0.1;
const double L = 1.479; 
const int N = int(L/DELTA)+1;
const deouble steps = 1000;


typedef std::vector<double> Matrix;

void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void print(const Matrix & m)
