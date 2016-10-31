#include <stdio.h>    // printf, fopen, fclose, fprintf, 
#include <stdlib.h>   // strtod 
#include <fstream>    // ifstream
#include <string>     // string
#include <string.h>
#include <ctime>      // clock

#define N 1048576 // 1024*1024 number of lines
//#define N 32

using namespace std;

void write(string name) {
  FILE* f = fopen(name.c_str(),"w");
  for(float i=0;i<N;i++) 
    fprintf(f,"  %s %4.2f %4.2f %4.2f\n","x",i,i,i); // write some data 
  fclose(f);
}

void read1(string name) {
  double num,check=0;
  string s;
  ifstream f(name);
  for(int i=0;i<N;i++) {
    f >> s;
    f >> num;
    f >> num;
    f >> num;
    check+=num;
  }
  printf("check = %.2f\n",check);
  f.close();
}

void read2(string name) {
  double num,check=0;
  char c[16];
  string s;
  FILE* f=fopen(name.c_str(),"r");
  while(fscanf(f,"%s%lf%lf%lf",c,&num,&num,&num)!=EOF) {
    s = c;
    check+=num;
  }
  printf("check = %.2f\n",check);
  fclose(f);
}

void read3(string name) {
  string line, s;
  double num,check=0;
  ifstream f(name);
  while(getline(f,line)) {
    size_t start = line.find_first_not_of(" \t");
    size_t pos = line.find(" ");
    char* c = &*(line.begin() + pos + 1);
    s = line.substr(start,pos+1); 
    num = strtod(c+start, &c);
    num = strtod(c, &c);
    num = strtod(c, &c);
    check+=num;
  }
  printf("check = %.2f\n",check);
  f.close();
}

int main() {
  clock_t start, end;
  string name("testfile.dat");

  start = clock();
  write(name);
  end   = clock();
  printf("write                 : took %.2f seconds\n",double(end-start)/CLOCKS_PER_SEC);

  start = clock();
  read1(name);
  end   = clock();
  printf("read1 (iostream)      : took %.2f seconds\n",double(end-start)/CLOCKS_PER_SEC);

  start = clock();
  read2(name);
  end   = clock();
  printf("read2 (fscanf)        : took %.2f seconds\n",double(end-start)/CLOCKS_PER_SEC);

  start = clock();
  read3(name);
  end   = clock();
  printf("read3 (stream+strtod) : took %.2f seconds\n",double(end-start)/CLOCKS_PER_SEC);
}
