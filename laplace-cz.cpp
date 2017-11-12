#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "papi.h"

//para compilar usar: 
// g++ -I $HOME/local/include -Wl,-rpath=$HOME/local/lib -std=c++11 laplace-cz.cpp -fopenmp -L $HOME/local/lib -lpapi


// constants
const double DELTA = 0.01;
const double L = 10.0; 
const int N = int(L/DELTA)+1;
const int STEPS = 100;
typedef std::vector<double> Matrix;

//Definición Funciones
void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void evolve(Matrix & m);
void print(const Matrix & m);
void init_gnuplot(void);
void plot_gnuplot(const Matrix & m);
void PAPI_measure_Init(int retval, float & ireal_time, float & iproc_time, float & imflops, long long & iflpops);
void PAPI_measure_Final(int retval, float & real_time, float & proc_time, float & mflops, long long & flpops);



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


  // ii=0;  //barra a portencial 100
  ii = N*0.25; // posicion vertical de barra
  for (jj = int(N*0.25); jj <int(N*0.75)+1; ++jj) //posicion horizontal
    m[ii*N + jj] = 100;

  ii = N*0.75; //barra a potencial -100 //posición Vertical
  for (jj = int(N*0.25); jj < int(N*0.75)+1; ++jj) //posicion horizontal
    m[ii*N + jj] = -100;

  
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


void evolve(Matrix & m)
{

  
#pragma omp parallel for
  for(int ii=0; ii<N; ++ii) {
    for(int jj=0; jj<N; ++jj) {
      // check if boundary
      if(ii == 0) continue;
      if(ii == int(N*0.25))continue; //barra a potencial 100
      if(ii == int(N*0.75))continue; // barra a potencial -100
      if(ii == N-1) continue;
      if(jj == 0) continue;
      if(jj == N-1) continue;
      // evolve non boundary
      m[ii*N+jj] = (m[(ii+1)*N + jj] +
                    m[(ii-1)*N + jj] +
                    m[ii*N + jj + 1] +
                    m[ii*N + jj - 1] )/4.0;
    }
  }
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


void init_gnuplot(void)
{
  std::cout << "set contour " << std::endl;
  std::cout << "set terminal gif animate " << std::endl;
  std::cout << "set out 'anim.gif' " << std::endl;
}

void plot_gnuplot(const Matrix & m)
{
  std::cout << "splot '-' w pm3d " << std::endl;
  print(m);
  std::cout << "e" << std::endl;  
}


void PAPI_measure_Init(int retval, float & ireal_time, float & iproc_time, float & imflops, long long & iflpops)
{
  if((retval=PAPI_flops(&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK)
    {
      printf("Could not initialise PAPI_flops \n");
      printf("Your platform may not support floating point operation event.\n");
      printf("retval: %d\n", retval);
      exit(1);
    }
}

void PAPI_measure_Final(int retval, float & real_time, float & proc_time, float & mflops, long long & flpops)
{
  if((retval=PAPI_flops( &real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
      printf("retval: %d\n", retval);
      exit(1);
    }
}
