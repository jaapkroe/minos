 awk 'BEGIN{
   pi2=2*atan2(0,-1); 
   N=50
   for(i=0;i<N;i++) {
     x=i/N; 
     # some perdiodic function f(i) such that f(0)=f(1)

     #print -1+cos(x*pi2)
     #print 1 + 2*sin(2*x*pi2) + 4*cos(4*x*pi2) 

     # or something that is not at all periodioc

     # like a Gaussian
     print exp(-(x-0.5)^2 / (2*0.1^2) ) 
     # or a polynomial
     #print 0.1 + 0.1 * (x-.4) - 2 * (x-.6)^2
     # some divergent function
     #print 1/(x-.230001)
   }
}' > fft_1d.dat
