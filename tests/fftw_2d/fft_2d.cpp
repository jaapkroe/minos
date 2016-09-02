#include <iostream>
#include <math.h>
#include <fstream>
#include <fftw3.h>
int main()
{
  int i, j, k, N1=0, N2=0, M1=128, M2=128;
  //M1=16, M2=16;
  char tmp[20];
  std::ifstream file("fft_2d.dat");
  file >> tmp >> N1 >> N2;
  // use complex numbers for input
  fftw_complex *input = fftw_alloc_complex(N1*N2);
  for(i=0;i<N1;i++) for(j=0;j<N2;j++) {
    int index = i*N2+j;
    file >> k >> k >> input[index][0];
    // set imaginary part to zero
    input[index][1]=0;
  }
  if(M1>N1)M1=N1; 
  if(M2>N2)M2=N2;
  //fprintf(stderr,"read %d x %d data points, generating FFT.\n",N1,N2);
  fprintf(stderr,"read %d x %d data points, generating FFT with %d x %d q-points\n",N1,N2,M1,M2);

  /* the following three lines are the core of the program */
  // allocate output, note that fft could also be done 'in-place' using the input twice
  fftw_complex *output = fftw_alloc_complex(N1*N2);
  // set up the fft, FFTW_FOWARD means the sign in exponent is positive
  // FFTW_ESIMATE is used for a basic fft, if the same plan is used many times FFTW_MEASURE would be a better choice
  fftw_plan p = fftw_plan_dft_2d(N1, N2, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  /*********************************************************/

  double norm = 1./(N1*N2);
  for(int qi=0;qi<N1;qi++) {
    for(int qj=0;qj<N2;qj++) {
      int index = qi*N2+qj;
      /* the follow line could be used instead to center the qi around zero */
      //printf("%-4d %-4d %12.4g %12.4g\n",(qi+N1/2)%N1-N1/2,(qj+N2/2)%N2-N2/2,output[index][0]*norm,output[index][1]*norm);
      printf("%-4d %-4d %12.4g %12.4g\n",qi,qj,output[index][0]*norm,output[index][1]*norm);
    }
    printf("\n");
  }

  fprintf(stderr,"generating backtransform to real data.\n");
  std::string fname="fft_2d.gen";
  FILE *fgen=fopen(fname.c_str(),"w");
  for(i=0;i<N1;i++) {
    double x = 2*M_PI*i/N1;
    for(j=0;j<N2;j++) {
      double y = 2*M_PI*j/N2;
      double f = 0;

      //for(int qi=0;qi<N1;qi++) {
      //  for(int qj=0;qj<N2;qj++) {
      
      /* the loops above are reformulated to pick a_q around q=0 
       * which is has periodic boundary conditions
       * they could be replaced by the commented-out loops in case 
       * one doesn't want to restrict the number of output basis functions 
       * for the expansion back.
       * The limits M1 and M2 serve only to illustrate how one can typically
       * recover rather well are 
       */
      for(int qi=-M1/2;qi<M1/2;qi++) {
				int qqi = (qi<0?qi+N1:qi);
        for(int qj=-M1/2;qj<M2/2;qj++) {
					int qqj = (qj<0?qj+N2:qj);
          int index = qqi*N2+qqj;

          /* needed if loops above are used instead ... */
          //int index = qi*N2+qj;
          
          f +=   output[index][0] * cos(qi*x + qj*y) \
               - output[index][1] * sin(qi*x + qj*y);
        }
      }
      fprintf(fgen,"%-4d %-4d %12.4g\n",i,j,f*norm);
    }
    fprintf(fgen,"\n");
  }
  fclose(fgen);

  /* for completeness we call the destructors explicity 
   * deallocate the memory and destroy the fftw plan
   */
  fftw_destroy_plan(p);
  fftw_free(input);
  fftw_free(output);
}
