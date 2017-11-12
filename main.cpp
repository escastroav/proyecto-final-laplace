#include "laplace.h"

int main(int argc, char ** argv)  //recive el número de nucleos a usar
{
  const int nucleos = std::atoi(argv[1]);
  Matrix data(N*N);

  float real_time, proc_time,mflops =0.0 ;
  long long flpops =0.0;
  float ireal_time , iproc_time, imflops =0.0 ;
  long long iflpops=0.0;
  int retval=0; 
 
  
  //Funciones para definir condiciones del espacio a estudiar  
  initial_conditions(data);
  boundary_conditions(data);
  //init_gnuplot();
  omp_set_num_threads(nucleos);

  for (int istep = 0; istep < STEPS; ++istep) {
    //LLAMADO DE FUNCIONES DENTRO DE PAPI
    PAPI_measure_Init(retval, ireal_time, iproc_time, imflops, iflpops);
    evolve(data); //funcion de evolución del sistema
    PAPI_measure_Final(retval, real_time, proc_time, mflops, flpops);   
    //plot_gnuplot(data);
    std::cout<<istep<<"\t"<<real_time<<"\t"<<proc_time<<"\t"<<flpops<<"\t"<<mflops<<std::endl;
    
  }
  
  return 0;
}
