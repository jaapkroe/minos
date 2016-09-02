 awk 'BEGIN{
   pi2=2*atan2(0,-1); 
   N1=64
   N2=64
   print "# ",N1,N2
   for(i=0;i<N1;i++) {
     x=i/N1; 
     for(j=0;j<N2;j++) {
       y=j/N2;
       #f = 1
       #f = cos(y*pi2)
       #f = cos(x*pi2) + 2*sin(2*y*pi2)
       #f = cos(x*pi2) + 2*sin(2*x*pi2) * cos(y*pi2)
       #f = cos(x*pi2)*cos(y*pi2)
       #f = 0.1*x**2-(y-.3)**3
       #f = exp( -((x-0.5)**2 *(2*y-0.3)**2 / (2*0.1**2)))
       f = exp( -((x-0.5)**2/(2*0.1**2) +(y-0.4)**2/(2*0.15**2)) )
       printf "%-4d %-4d %12.4g\n",i,j,f
     }
     print ""
   }
}' > fft_2d.dat
