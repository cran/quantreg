C Output from Public domain Ratfor, version 1.05
      subroutine pfnb(n,p,m,a,y,q,r,b,band,m0,d,u,wn,wp, aa,yy,slo,shi,r
     *hs,glob,ghib,nit,info)
      integer n,p,m,kk(2),mm,m0,nit(5,m),info,sumbad
      integer i,j,k,slo(n),shi(n),ifix,ibad
      logical notopt
      double precision a(p,n),y(n),r(n),q(m),d(n),u(n),b(p,m)
      double precision wn(n,9), wp(p,(p+3)),band(n),qk(2)
      double precision glob(p),ghib(p),aa(p,n),yy(n),rhs(p)
      double precision zero,one,tau,beta,eps,big,fm,fn
      parameter(zero = 0.0d0)
      parameter(one = 1.0d0)
      parameter(beta = 0.99995d0)
      parameter(big = 1.0d+10)
      parameter(eps = 1.0d-06)
      do23000 iq = 1,m 
      notopt = .true.
      tau = q(iq)
      mm = m0
      ifix = 0
      ibad = 0
23002 if(notopt)then
      ibad = ibad + 1
      fm = mm
      fn = dble(n)
      kk(1) = int(n * dmax1(1./fn, tau - fm/(2 * fn))) + 1
      kk(2) = int(n * dmin1(tau + fm/(2. * fn), (fn - 1)/fn))
      do23004 i = 1,n
      u(i) = r(i)/band(i)
23004 continue
23005 continue
      call kuantiles(kk,2,n,u)
      qk(1) = u(kk(1))
      qk(2) = u(kk(2))
      call iphil(n,0,slo)
      call iphil(n,0,shi)
      do23006 i = 1,n
      if(r(i) .lt. (band(i) * qk(1)))then
      slo(i) = 1
      else
      if(r(i) .gt. (band(i) * qk(2)))then
      shi(i) = 1
      endif
      endif
23006 continue
23007 continue
23012 if(notopt)then
      ifix = ifix + 1
      call dphil(p,zero,glob)
      call dphil(p,zero,ghib)
      call dphil(n,one,d)
      call dphil(n,one,u)
      k = 0
      do23014 i = 1,n
      if(slo(i) .eq. 0 .and. shi(i) .eq. 0)then
      k = k + 1
      call dcopy(p,a(1,i),1,aa(1,k),1)
      yy(k) = -y(i)
      else
      if(slo(i) .eq. 1)then
      do23020 j = 1,p
      glob(j) = glob(j) + a(j,i)
23020 continue
23021 continue
      else
      if(shi(i) .eq. 1)then
      do23024 j = 1,p
      ghib(j) = ghib(j) + a(j,i)
23024 continue
23025 continue
      endif
      endif
      endif
23014 continue
23015 continue
      call dcopy(p,glob,1,aa(1,k+1),1)
      call dcopy(p,ghib,1,aa(1,k+2),1)
      yy(k+1) = big
      yy(k+2) = -big
      call dgemv('N',p,k+2,one-tau,aa,p,d,1,zero,rhs,1)
      call dscal(k+2,zero,wn,1)
      call daxpy(k+2,one-tau,u,1,wn,1)
      call rqfnb(k+2,p,aa,yy,rhs,d,u,beta,eps,wn,wp,nit(1,iq),info)
      call dcopy(p,wp,1,b(1,iq),1)
      call dcopy(n,y,1,r,1)
      call dgemv('T',p,n,one,a,p,b(1,iq),1,one,r,1)
      sumbad = 0
      do23026 i = 1,n 
      if((r(i) .gt. 0) .and. slo(i) .eq. 1)then
      slo(i) = 0
      sumbad = sumbad + 1
      endif
      if((r(i) .lt. 0) .and. shi(i) .eq. 1)then
      shi(i) = 0
      sumbad = sumbad + 1
      endif
23026 continue
23027 continue
      if(sumbad .gt. 0)then
      if(sumbad .gt. 0.1 * mm)then
      mm = min(2 * mm, n)
      goto 23013
      endif
      else
      notopt = .false.
      endif
      goto 23012
      endif
23013 continue
      nit(4,iq) = ifix
      nit(5,iq) = ibad
      goto 23002
      endif
23003 continue
23000 continue
23001 continue
      return
      end
      subroutine iphil(n,a,v)
      integer n,i,a,v(n)
      do23036 i = 1,n
      v(i) = a
23036 continue
23037 continue
      return
      end
      subroutine dphil(n,a,v)
      integer i,n
      double precision a,v(n)
      do23038 i = 1,n
      v(i) = a
23038 continue
23039 continue
      return
      end
