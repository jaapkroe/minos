#include <stdio.h>
#include <random>
#include <getopt.h>
#include <cmath>
#include <string>
using namespace std;

int main(int argc, char** argv) {
  char c;
  int n=1000;
  double rho=0.2;
  bool random=false;
  string usage="usage: write-config -n <natoms> -d <density> [-r for randomized coordinates] [-h for help]\n";

  while ((c = getopt(argc, argv, "d:n:rh")) != -1) {
    switch(c) {
      case 'n': n=(int)atof(optarg); break;
      case 'd': rho=atof(optarg); break;
      case 'r': random=!random; break;
      case 'h': fprintf(stderr,"%s",usage.c_str()); return 0;
      default : fprintf(stderr,"%s",usage.c_str()); return 1;
    }
  }

  FILE *f;
  f = fopen("config.xyz","w");
  // n = rho * l^3 -> l = (n/rho)^1/3
  if(random) {
    double l=pow(float(n)/float(rho),1./3.);
    fprintf(f,"%d\n%f %f %f\n",n,l,l,l);
    for(int i=0;i<n;i++) 
      fprintf(f,"C %.2f %.2f %.2f\n",rand()*l/RAND_MAX,rand()*l/RAND_MAX,rand()*l/RAND_MAX);
  } else {
    int ni = pow(n,1/3.)+1;
    double l=float(ni)/float(rho)/3.;
    double d=l/ni;
    fprintf(f,"%d %d\n%f %f %f\n",n,ni,l,l,l);
    int counter=0;
    for(int i=0;i<ni;i++) {
      for(int j=0;j<ni;j++) {
        for(int k=0;k<ni;k++,counter++) {
          if(counter>=n) continue;
          fprintf (f,"C %.2f %.2f %.2f\n",i*d,j*d,k*d);
        }
      }
    }
  }

  fclose(f);
}

