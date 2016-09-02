#include <iostream>
#include <fstream>
#include <fftw3.h>
#define NMAX 8192
int main()
{
  int i, N=0, M=10;
  double *input = new double [NMAX];
  std::ifstream file("fft_1d.dat");
  while (file >> input[N]) { N++; if(N>=NMAX) { fprintf(stderr,"ERROR: N>NMAX\n"); exit(1); }; };
  if(M>N) M=N;
  fprintf(stderr,"read %d data points, generating FFT with %d q-points\n",N,M);

  fftw_complex *out = fftw_alloc_complex(2*N/1);
  fftw_plan p = fftw_plan_dft_r2c_1d(N,input,out,FFTW_ESTIMATE);
  fftw_execute(p);

  for(i=0;i<M;i++)
    printf("%-4d %12.4g %12.4g\n",i,out[i][0]/N,out[i][1]/N);
}
