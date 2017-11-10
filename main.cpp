#include "laplace.h"

int main(void)
{
  Matrix data(N*N);
  initial_conditions(data);
  boundary_conditions(data);

  //evolve(data);
  //print(data);
  
  init_gnuplot();
  for (int istep = 0; istep < STEPS; ++istep) {
    evolve(data);
    plot_gnuplot(data);
  }
  
  return 0;
}
