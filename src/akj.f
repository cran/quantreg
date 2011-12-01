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
     *     one,two
      parameter( zero = 0.d0)
      parameter( one = 1.d0)
      parameter( two = 2.d0)
      parameter( four = 4.d0)
      parameter( half = 0.5d0)
      parameter( fifth = 0.2d0)
      parameter( pi = 3.141593d0)

      xn=nx

C call srtad(x,1,nx) #port sort routine now done in S interface.

      if(iker.eq.0) then
         con1= one/sqrt(2.0*pi)
      else if(iker.eq.1) then
         con1= one/pi
      endif

C if no h is provided, calculate a default
      if(h.le.0.) then
         sum=0.
         sqsum=0.
         do 23006 i=1,nx
            sqsum=sqsum+x(i)*x(i)*p(i)
            sum=sum+x(i)*p(i)
23006    continue
         xsd=dsqrt(sqsum-sum*sum)
c     compute 'qrange' :=  IQR (x[i])

         sum=zero
c        first, qrange = Q_1 [ = quantile(*, 0.25) ]
         do i=1,nx
            sum=sum+p(i)
            if(sum .ge. .25) then
               qrange = x(i)
               goto 23010
            endif
         enddo
23010    continue
         sum=one
         do i=nx,1,-1
            sum=sum-p(i)
            if(sum .le. .75) then
               qrange = x(i) - qrange
               goto 23015
            endif
         enddo
23015    continue

         a=min(xsd,qrange/1.34)
         h=kappa*a/(xn**fifth)
C                             see Silverman p 48
      endif
      hinv=one/h

C Stage one:  compute pilot estimate of density
      do 23018 j=1,nx
         xker=0.
         if(iker.eq.0) then
            do 23022 i=1,nx
               xponen=(x(j)-x(i))*hinv
               xponen=half*xponen**2
               xker=xker+p(i)*exp(-xponen)*hinv
23022       continue
         else if(iker.eq.1) then
            do 23026 i=1,nx
               xponen=(x(j)-x(i))*hinv
               xker=xker+p(i)*hinv/(1+xponen**2)
23026       continue
         endif
         xlam(j)=con1*xker
23018 continue

C Stage two:  Automatic window widths (Silverman p101)
      glog=zero
      do 23028 i=1,nx
         glog=glog+p(i)*log(xlam(i))
23028 continue
      g=exp(glog)
      ginv=one/g
      do 23030 i=1,nx
         xlam(i)=hinv/((xlam(i)*ginv)**(-alpha))
C notice xlam no longer its own self at this pt! xlam is 1/(h*lambda(i))
C substitution of * for / thus achieved speeds things up
C Stage two:  new density-score estimates
23030 continue

      do 23032 j=1,nz
         xker=zero
         dxker=zero
         ddxker=zero
         if(iker.eq.0) then
C  gaussian kernel

            do 23036 i=1,nx
               xponen=(z(j)-x(i))*xlam(i)
               fact=exp(-half*xponen*xponen)*xlam(i)
               xker=xker+p(i)*fact
               dxker=dxker-p(i)*fact*xponen*xlam(i)
               ddxker=ddxker- p(i)*fact*(one - xponen**2)*xlam(i)**2
23036       continue

         else if(iker.eq.1) then
C   cauchy kernel

            do 23040 i=1,nx
               xponen=(z(j)-x(i))*xlam(i)
               fact=xlam(i)/(one+xponen**2)
               xker=xker+p(i)*fact
               dxker=dxker-p(i)*two*xponen*fact**2
               ddxker=ddxker- p(i)*two*(fact**2)*
     *              (xlam(i)- four*(xponen**2)*fact)
23040       continue
         endif
         dens(j)=con1*xker
         psi(j)=-(dxker/xker)
         score(j)=(dxker/xker)**2-ddxker/xker
23032 continue

      return
      end
