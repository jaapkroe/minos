#include <stdio.h>
#include <random>
#include <getopt.h>
#include <cmath>
#include <string>
using namespace std;

int main(int argc, char** argv) {
  char c;
  int n=1000;
  double rho=1;
  string usage="usage: write-config -n <natoms> -d <density>\n";

  while ((c = getopt(argc, argv, "d:n:h")) != -1) {
    switch(c) {
      case 'n': n=atoi(optarg); break;
      case 'd': rho=atof(optarg); break;
      case 'h': fprintf(stderr,"%s",usage.c_str()); return 0;
      default : fprintf(stderr,"%s",usage.c_str()); return 1;
    }
  }

  FILE *f;
  f = fopen("config.xyz","w");
  // n = rho * l^3 -> l = (n/rho)^1/3
  double l=pow(float(n)/float(rho),1./3.);
  fprintf(f,"%d\n%f %f %f\n",n,l,l,l);
  for(int i=0;i<n;i++) 
    fprintf(f,"C %.2f %.2f %.2f\n",rand()*l/RAND_MAX,rand()*l/RAND_MAX,rand()*l/RAND_MAX);
  fclose(f);
}

