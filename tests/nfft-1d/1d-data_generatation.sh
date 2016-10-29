 awk 'BEGIN{
   pi2=2*atan2(0,-1); 
   N=16
   M=64
   L=1.0
   print "# ",M,N,L
   for(i=0;i<M;i++) {
     x=i/M; 
     # some perdiodic function f(i) such that f(0)=f(1)

     f=10*sin(x*pi2)
     #f=-1+cos(x*pi2)
     #f=1 + 2*sin(2*x*pi2) + 4*cos(4*x*pi2) 

     # or something that is not at all periodioc

     # like a Gaussian
     #f=exp(-(x-0.5)^2 / (2*0.1^2) ) 
     # or a polynomial
     #f=0.1 + 0.1 * (x-.4) - 2 * (x-.6)^2
     # some divergent function
     #f=1/(x-.230001)

     # imaginary part
     fimg=0.0
     printf("%10.4g %12.6g %12.6g\n",x,f,fimg);
   }
}' > nfft_1d.dat
