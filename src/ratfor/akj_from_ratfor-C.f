C Output from Public domain Ratfor, version 1.01
      subroutine akj(x,z,p,iker,dens,psi,score,nx,nz,h,alpha,kappa,xlam)
C univariate kernel density-score estimator
C the algorithm is basically from Silverman as adapted for Portnoy and Koenker
C Annals paper on adaptive L-estimation in regression.
C x--pts used for centers of kernel assumed to be sorted!!!!
C z--pts at which density is calculated
C p--probability associated with x's
C dens--f(z), the densities at z
C psi--f'(z)/f(z) the score at z
C score--(log(f(z)))'', the J fn at z
C nx--number of pts in x
C nz--number of pts in z
C iker--kernel 
C          0=gaussian
C	   1=cauchy
C h--initial window size (overall)--choose zero for default
C kappa--constant determining initial (default) window width
C xlam--Silverman's lambda, window adjustment for each x
      double precision dens(nz),score(nz),psi(nz),h,kappa
      double precision z(nz),x(nx),xlam(nx),p(nx),qrange,pi
      double precision con1,sum,sqsum,xsd,a,fifth,hinv,half
      double precision xn,xker,dxker,ddxker,fact,xponen,alpha,glog,zero,
     *one,two
      parameter( zero = 0.d0)
      parameter( one = 1.d0)
      parameter( two = 2.d0)
      parameter( four = 4.d0)
      parameter( half = 0.5d0)
      parameter( fifth = 0.2d0)
      parameter( pi = 3.141593d0)
      xn=nx
Ccall srtad(x,1,nx) #port sort routine now done in S interface.
      if(iker.eq.0)then
      con1=one/sqrt(2.0*pi)
      else
      if(iker.eq.1)then
      con1=one/pi
Cif no h is provided, calculate a default
      endif
      endif
      if(h.le.0.)then
      sum=0.
      sqsum=0.
      do23006 i=1,nx
      sqsum=sqsum+x(i)*x(i)*p(i)
      sum=sum+x(i)*p(i)
23006 continue
23007 continue
      xsd=dsqrt(sqsum-sum*sum)
      sum=zero
      i=1
23008 if(.not.(i.lt.nx))goto 23010
      sum=sum+p(i)
      if(sum.lt..25)then
      goto 23009
      else
      qrange=x(i)
      goto 23010
      endif
23009 i=i+1
      goto 23008
23010 continue
      sum=one
      i=nx
23013 if(.not.(i.gt.0))goto 23015
      sum=sum-p(i)
      if(sum.gt..75)then
      goto 23014
      else
      qrange=x(i)-qrange
      goto 23015
      endif
23014 i=i-1
      goto 23013
23015 continue
      a=min(xsd,qrange/1.34)
Csee Silverman p 48
      h=kappa*a/(xn**fifth)
      endif
      hinv=one/h
CStage one:  compute pilot estimate of density
      do23018 j=1,nx
      xker=0.
      if(iker.eq.0)then
      do23022 i=1,nx 
      xponen=(x(j)-x(i))*hinv
      xponen=half*xponen**2
      xker=xker+p(i)*exp(-xponen)*hinv
23022 continue
23023 continue
      else
      if(iker.eq.1)then
      do23026 i=1,nx 
      xponen=(x(j)-x(i))*hinv
      xker=xker+p(i)*hinv/(1+xponen**2)
23026 continue
23027 continue
      endif
      endif
      xlam(j)=con1*xker
CStage two:  Automatic window widths (Silverman p101)
23018 continue
23019 continue
      glog=zero
      do23028 i=1,nx
      glog=glog+p(i)*log(xlam(i))
23028 continue
23029 continue
      g=exp(glog)
      ginv=one/g
      do23030 i=1,nx
      xlam(i)=hinv/((xlam(i)*ginv)**(-alpha))
Cnotice xlam no longer its own self at this pt! xlam is 1/(h*lambda(i))
Csubstitution of * for / thus achieved speeds things up
CStage two:  new density-score estimates
23030 continue
23031 continue
      do23032 j=1,nz
      xker=zero
      dxker=zero
      ddxker=zero
      if(iker.eq.0)then
      do23036 i=1,nx 
Cgaussian kernel
      xponen=(z(j)-x(i))*xlam(i)
      fact=exp(-half*xponen*xponen)*xlam(i)
      xker=xker+p(i)*fact
      dxker=dxker-p(i)*fact*xponen*xlam(i)
      ddxker=ddxker- p(i)*fact*(one - xponen**2)*xlam(i)**2
23036 continue
23037 continue
      endif
      if(iker.eq.1)then
C cauchy kernel 
      do23040 i=1,nx 
      xponen=(z(j)-x(i))*xlam(i)
      fact=xlam(i)/(one+xponen**2)
      xker=xker+p(i)*fact
      dxker=dxker-p(i)*two*xponen*fact**2
      ddxker=ddxker- p(i)*two*(fact**2)*(xlam(i)- four*(xponen**2)*fact)
23040 continue
23041 continue
      endif
      dens(j)=con1*xker
      psi(j)=-(dxker/xker)
      score(j)=(dxker/xker)**2-ddxker/xker
23032 continue
23033 continue
      return
      end
