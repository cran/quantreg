      subroutine akj(x,z,p,iker,dens,psi,score,nx,nz,h,alpha,kappa,xlam)
      double precision dens(nz),score(nz),psi(nz),h,kappa
      double precision z(nz),x(nx),xlam(nx),p(nx),qrange,pi
      double precision con1,sum,sqsum,xsd,a,fifth,hinv,half
      double precision xn,xker,dxker,ddxker,fact,xponen,alpha,glog,zero,
&     one,two
      data zero/0.d0/
      data one/1.d0/
      data two/2.d0/
      data four/4.d0/
      data half/0.5d0/
      data fifth/0.2d0/
      data pi/3.141593d0/
      xn=nx
      if(.not.(iker.eq.0))goto 23000
      con1=one/sqrt(2.0*pi)
      goto 23001
23000 continue
      if(.not.(iker.eq.1))goto 23002
      con1=one/pi
23002 continue
23001 continue
      if(.not.(h.le.0.))goto 23004
      sum=0.
      sqsum=0.
      do 23006 i=1,nx
      sqsum=sqsum+x(i)*x(i)*p(i)
      sum=sum+x(i)*p(i)
23006 continue
      xsd=dsqrt(sqsum-sum*sum)
      sum=zero
      i=1
23008 if(.not.(i.lt.nx))goto 23010
      sum=sum+p(i)
      if(.not.(sum.lt..25))goto 23011
      goto 23009
23011 continue
      qrange=x(i)
      goto 23010
23012 continue
23009 i=i+1
      goto 23008
23010 continue
      sum=one
      i=nx
23013 if(.not.(i.gt.0))goto 23015
      sum=sum-p(i)
      if(.not.(sum.gt..75))goto 23016
      goto 23014
23016 continue
      qrange=x(i)-qrange
      goto 23015
23017 continue
23014 i=i-1
      goto 23013
23015 continue
      a=min(xsd,qrange/1.34)
      h=kappa*a/(xn**fifth)
23004 continue
      hinv=one/h
      do 23018 j=1,nx
      xker=0.
      if(.not.(iker.eq.0))goto 23020
      do 23022 i=1,nx 
      xponen=(x(j)-x(i))*hinv
      xponen=half*xponen**2
      xker=xker+p(i)*exp(-xponen)*hinv
23022 continue
      goto 23021
23020 continue
      if(.not.(iker.eq.1))goto 23024
      do 23026 i=1,nx 
      xponen=(x(j)-x(i))*hinv
      xker=xker+p(i)*hinv/(1+xponen**2)
23026 continue
23024 continue
23021 continue
      xlam(j)=con1*xker
23018 continue
      glog=zero
      do 23028 i=1,nx
      glog=glog+p(i)*log(xlam(i))
23028 continue
      g=exp(glog)
      ginv=one/g
      do 23030 i=1,nx
      xlam(i)=hinv/((xlam(i)*ginv)**(-alpha))
23030 continue
      do 23032 j=1,nz
      xker=zero
      dxker=zero
      ddxker=zero
      if(.not.(iker.eq.0))goto 23034
      do 23036 i=1,nx 
      xponen=(z(j)-x(i))*xlam(i)
      fact=exp(-half*xponen*xponen)*xlam(i)
      xker=xker+p(i)*fact
      dxker=dxker-p(i)*fact*xponen*xlam(i)
      ddxker=ddxker-p(i)*fact*(one - xponen**2)*xlam(i)**2
23036 continue
23034 continue
      if(.not.(iker.eq.1))goto 23038
      do 23040 i=1,nx 
      xponen=(z(j)-x(i))*xlam(i)
      fact=xlam(i)/(one+xponen**2)
      xker=xker+p(i)*fact
      dxker=dxker-p(i)*two*xponen*fact**2
      ddxker=ddxker-p(i)*two*(fact**2)*(xlam(i)-four*(xponen**2)*fact)
23040 continue
23038 continue
      dens(j)=con1*xker
      psi(j)=-(dxker/xker)
      score(j)=(dxker/xker)**2-ddxker/xker
23032 continue
      return
      end
