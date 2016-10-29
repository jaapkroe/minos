#include <fftw3.h>
#include <nfft3.h>

//  void simple_test_nfft_1d()
//  {
//    int j,k;
//    nfft_plan my_plan;
//  
//    nfft_init_1d(&my_plan,12,19);
//    
//    for(j=0;j<my_plan.M;j++)
//      my_plan.x[j]=((double)rand())/RAND_MAX-0.5;
//  
//    if(my_plan.nfft_flags & PRE_PSI)
//      nfft_precompute_psi(&my_plan);
//  
//    for(k=0;k<my_plan.N_L;k++)
//      {
//        my_plan.f_hat[k][0]=((double)rand())/RAND_MAX;
//        my_plan.f_hat[k][1]=((double)rand())/RAND_MAX;
//      }
//    
//    nfft_trafo(&my_plan);
//  
//    nfft_finalize(&my_plan);
//  }
//  

void simple_test_nfft_1d()
{
  nfft_plan p;
  int N=14;
  int M=19;

  nfft_init_1d(&p,N,M);

  nfft_vrand_shifted_unit_double(p.x,p.M_total);

  if(p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);

  nfft_vrand_unit_complex(p.f_hat,p.N_total);
  nfft_vpr_complex(p.f_hat,p.N_total,"given Fourier coefficients, f_hat");

  ndft_trafo(&p);
  nfft_vpr_complex(p.f,p.M_total,"ndft, f");

  nfft_trafo(&p);
  nfft_vpr_complex(p.f,p.M_total,"nfft, f");

  nfft_finalize(&p);
}
