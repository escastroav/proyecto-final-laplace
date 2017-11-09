//#include "laplace.h"
#include <iostream>
#include <vector>
#include <cmath>


// constants
const double DELTA = 0.1;
const double L = 1.479; 
const int N = int(L/DELTA)+1;
const double steps = 1000;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void print(const Matrix & m);

int main(void)
{
  Matrix data(N*N);
  initial_conditions(data);
  boundary_conditions(data);

  // evolve(data);
  
  print(data);
  
  return 0;
}

void initial_conditions(Matrix & m)
{
  for(int ii=0; ii<N; ++ii) {
    for(int jj=0; jj<N; ++jj) {
      m[ii*N + jj] = 1.0;
    }
  }
}

void boundary_conditions(Matrix & m)
{
  int ii = 0, jj = 0;

  ii = 0;
  for (jj = 0; jj < N; ++jj)
    m[ii*N + jj] = 100;

  ii = N-1;
  for (jj = 0; jj < N; ++jj)
    m[ii*N + jj] = 0;

  jj = 0;
  for (ii = 1; ii < N-1; ++ii)
    m[ii*N + jj] = 0;

  jj = N-1;
  for (ii = 1; ii < N-1; ++ii)
    m[ii*N + jj] = 0;
}

void print(const Matrix & m)
{
  for(int ii=0; ii<N; ++ii) {
    for(int jj=0; jj<N; ++jj) {
      std::cout << ii*DELTA << " " << jj*DELTA << " " <<  m[ii*N + jj] << "\n";
    }
    std::cout << "\n";
  }  
}
