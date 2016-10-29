#include "stdio.h"
#include "math.h"

void pack(unsigned long long &hash, unsigned i, unsigned j, unsigned k, unsigned n) {
  unsigned long long ii(i),jj(j),kk(k);
  hash = ( ii ) + ( jj << (n) ) + ( kk << (2*n) );  
  // bitwise shift :
  // counting from right, the first n bits are for i, then n for j, then n for k
  // so if this hash is 64-bit, n can be at most 21
}

void unpack(unsigned long long hash, unsigned &i, unsigned &j, unsigned &k, unsigned n) {
  unsigned long long len = pow(2,n)-1;
  //unsigned long long len = 65535;
  k = ( hash >> (2*n) ) & len;  // shift right by 2n, then take the last n digits in binary
  j = ( hash >> n ) & len;      // shift right by  n, then take the last n digits in binary
  i = hash & len;               // i is simply the last n binary digits
}

int main() {
  unsigned i,j,k;
  unsigned long long hash;
  unsigned intlen = 21;
  unsigned long long maxint = pow(2,intlen);
  int nbits = sizeof(hash);
  if(3*intlen>nbits*8) printf("WARNING : packing factor (3x%d) exceeds hash dimensions (%d) in bit-size. This will result in unreliable bit-packing and failures for large k.\n",intlen,nbits*8);

  i=21238;
  j=13;
  k=121;
  if(i>maxint) printf("WARNING : i exceeds the integer length for packing (%d)\n",maxint);
  if(j>maxint) printf("WARNING : j exceeds the integer length for packing (%d)\n",maxint);
  if(k>maxint) printf("WARNING : k exceeds the integer length for packing (%d)\n",maxint);

  printf("prepacking : i,j,k = %d,%d,%d\n",i,j,k);

  pack(hash,i,j,k, intlen);
  printf("             -> computed hash  = %lld\n",hash);

  unpack(hash,i,j,k, intlen);
  printf("unpacking  : i,j,k = %d,%d,%d\n",i,j,k);

}
