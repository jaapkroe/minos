 awk 'BEGIN{
   pi2=2*atan2(0,-1); 
   N=64
   L=1.0
   print "# ",N,L
   for(i=0;i<N;i++) {
     x=i/N; 
     # some perdiodic function f(i) such that f(0)=f(1)

     #f=-1+cos(x*pi2)
     #f=1 + 2*sin(2*x*pi2) + 4*cos(4*x*pi2) 

     # or something that is not at all periodioc

     # like a Gaussian
     f=exp(-(x-0.5)^2 / (2*0.1^2) ) 
     # or a polynomial
     #f= 0.1 + 0.1 * (x-.4) - 2 * (x-.6)^2
     # some divergent function
     #f= 1/(x-.230001)

     printf "%12.4f %16.8f\n",x*pi2,f
   }
}' > nfft_1d.dat
