/*
 * Copyright (c) 2002, 2016 Jens Keiner, Stefan Kunis, Daniel Potts
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define NFFT_PRECISION_DOUBLE

#include "nfft3mp.h"

static void simple_test_nfft_1d(void)
{
  NFFT(plan) p;

  int N = 14; // number of fourier coefficients
  int M = 24; // number of data points
  int i;
  double L;
  FILE* f;
  char c[16];
  const char *error_str;

  // read data file
  f = fopen("nfft_1d.dat","r");
  fscanf(f, "%s %d %d %lf", c, &M, &N, &L);                                   // number of points, coefficients and box size
  NFFT(init_1d)(&p, N, M);                                                    // initialize the FFT plan (p) and allocate arrays for p
  for(i=0;i<M;i++) fscanf(f, "%lf %lf %lf", &p.x[i], &p.f[i][0], &p.f[i][1]); // read data file
  fclose(f);
  for(i=0;i<M;i++) { p.x[i]/=L; p.x[i] -= round(p.x[i]); }                    // wrap data to cell [-0.5, 0.5)
  f=fopen("input-check.dat","w"); for(i=0;i<M;i++) fprintf(f,"%-3d %12.8g %12.8g %12.8g\n",i,p.x[i],p.f[i][0],p.f[i][1]); fclose(f);

  //  /** init pseudo random nodes */
  NFFT(vrand_shifted_unit_double)(p.x, p.M_total);

  /** precompute psi, the entries of the matrix B */
  if (p.flags & PRE_ONE_PSI) NFFT(precompute_one_psi)(&p);

  /** init pseudo random Fourier coefficients and show them */
  NFFT(vrand_unit_complex)(p.f_hat,p.N_total);
  NFFT(vpr_complex)(p.f_hat, p.N_total, "given Fourier coefficients, vector f_hat");

  /** check for valid parameters before calling any trafo/adjoint method */
  error_str = NFFT(check)(&p);
  if (error_str != 0)
  {
    printf("Error in nfft module: %s\n", error_str);
    return;
  }

  /** direct trafo and show the result */
  NFFT(trafo_direct)(&p);
  NFFT(vpr_complex)(p.f,p.M_total,"ndft, vector f");
  f=fopen("output_direct.dat","w"); for(i=0;i<M;i++) fprintf(f,"%-3d %12.8g %12.8g %12.8g\n",i,p.x[i],p.f[i][0],p.f[i][1]); fclose(f);

  // /** approx. trafo and show the result */
  NFFT(trafo)(&p);
  NFFT(vpr_complex)(p.f,p.M_total,"nfft, vector f");
  f=fopen("output_approx.dat","w"); for(i=0;i<M;i++) fprintf(f,"%-3d %12.8g %12.8g %12.8g\n",i,p.x[i],p.f[i][0],p.f[i][1]); fclose(f);

  // /** approx. adjoint and show the result */
  NFFT(adjoint_direct)(&p);
  NFFT(vpr_complex)(p.f_hat,p.N_total,"adjoint ndft, vector f_hat");
  f=fopen("output_adj_direct.dat","w"); for(i=0;i<N;i++) fprintf(f,"%-3d %12.8g %12.8g\n",i,p.f_hat[i][0],p.f_hat[i][1]); fclose(f);

  // /** approx. adjoint and show the result */
  NFFT(adjoint)(&p);
  NFFT(vpr_complex)(p.f_hat,p.N_total,"adjoint nfft, vector f_hat");
  f=fopen("output_adj_approx.dat","w"); for(i=0;i<N;i++) fprintf(f,"%-3d %12.8g %12.8g\n",i,p.f_hat[i][0],p.f_hat[i][1]); fclose(f);

  /** finalise the one dimensional plan */
  NFFT(finalize)(&p);
}

int main(void)
{
  printf("1) computing a one dimensional ndft, nfft and an adjoint nfft\n\n");
  simple_test_nfft_1d();

  return EXIT_SUCCESS;
}
