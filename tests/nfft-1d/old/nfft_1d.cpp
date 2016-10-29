#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <nfft3.h>
#define NMAX 8192
int main()
{
  int i, N, M=16;
  double L;
  char tmp[16];
  std::ifstream file("nfft_1d.dat");
  file >> tmp >> N >> L;
  double *input = new double[N];
  double *grid = new double[N];
  if(N>NMAX) { fprintf(stderr, "N>NMAX (%d>%d)\n",N,NMAX); exit(1); };
  for(i=0;i<N;i++) file >> grid[i] >> input[i];

  nfft_plan p;
  fprintf(stdout,"init_1d:\n");
  nfft_init_1d(&p,N,M);
  if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  fprintf(stdout,"vrand shifted unit double:\n");
  nfft_vrand_shifted_unit_double(p.x,p.M_total);
  fprintf(stdout,"vpr_complex:\n");
  nfft_vpr_complex(p.f_hat,p.N_total,"given Fourier coefficients, f_hat");

  fprintf(stdout,"trafo:\n");
  nfft_trafo(&p);

  fprintf(stdout,"finalize:\n");
  nfft_finalize(&p);

  //if(M>N) M=N;
  fprintf(stdout,"read %d data points, generating FFT with %d q-points\n",N,M);

  //fftw_complex *out = fftw_alloc_complex(M);
  //fftw_plan p = fftw_plan_dft_r2c_1d(N, input, out, FFTW_ESTIMATE);
  //fftw_execute(p);

  //for(i=0;i<M;i++) printf("%-4d %12.4g %12.4g\n",i,out[i][0]/N,out[i][1]/N);
}
