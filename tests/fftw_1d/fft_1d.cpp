#include <iostream>
#include <fstream>
#include <fftw3.h>
#define NMAX 8192
int main()
{
  int i, N, M=10;
  char tmp[16];
  std::ifstream file("fft_1d.dat");
  file >> tmp >> N;
  double *input = new double[N];
  for(i=0;i<N;i++) file >> input[i];
  if(M>N) M=N;
  fprintf(stderr,"read %d data points, generating FFT with %d q-points\n",N,M);

  fftw_complex *out = fftw_alloc_complex(N);
  fftw_plan p = fftw_plan_dft_r2c_1d(N, input, out, FFTW_ESTIMATE);
  fftw_execute(p);

  for(i=0;i<M;i++) printf("%-4d %12.4g %12.4g\n",i,out[i][0]/N,out[i][1]/N);
}
