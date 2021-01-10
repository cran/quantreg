C Output from Public domain Ratfor, version 1.05
      subroutine pwxy(n,p,m,a,y,tau,qk,r,b,w,band,n0,d,u,wn,wp, aa,yy,sl
     *o,shi,rhs,glob,ghib,nit,info)
      integer n,p,m,kk(2),nn,n0,nit(5,m),info(m),sumbad,i,j,k,ir
      integer loq,hiq, slo(n),shi(n),ifix,ibad
      logical notopt
      double precision a(p,n),y(n),tau,qk(2),r(n),d(n),u(n),b(p,m),w(n)
      double precision wn(n,9), wp(p,(p+3)),band(n)
      double precision glob(p),ghib(p),aa(p,n),yy(n),rhs(p)
      double precision zero,one,beta,eps,big
      parameter(zero = 0.0d0)
      parameter(one = 1.0d0)
      parameter(beta = 0.99995d0)
      parameter(big = 1.0d+10)
      parameter(eps = 1.0d-06)
      do23000 ir = 1,m 
      notopt = .true.
      call grexp(n,w,one)
      nn = n0
      ifix = 0
      ibad = 0
23002 if(notopt)then
      ibad = ibad + 1
      loq = max0(1, int(n*tau - nn/2.))
      hiq = min0(int(n*tau + nn/2.), n)
      qk(1) = r(loq)
      qk(2) = r(hiq)
      call iphil(n,0,slo)
      call iphil(n,0,shi)
      do23004 i = 1,n
      if(r(i) .lt. qk(1))then
      slo(i) = 1
      else
      if(r(i) .gt. qk(2))then
      shi(i) = 1
      endif
      endif
23004 continue
23005 continue
23010 if(notopt)then
      ifix = ifix + 1
      call dphil(p,zero,glob)
      call dphil(p,zero,ghib)
      call dphil(n,one,d)
      call dphil(n,one,u)
      k = 0
      do23012 i = 1,n
      if(slo(i) .eq. 0 .and. shi(i) .eq. 0)then
      k = k + 1
      call dphil(p,zero,aa(1,k))
      call daxpy(p,w(i),a(1,i),1,aa(1,k),1)
      yy(k) = -y(i)*w(i)
      else
      if(slo(i) .eq. 1)then
      do23018 j = 1,p
      glob(j) = glob(j) + a(j,i) * w(i)
23018 continue
23019 continue
      else
      if(shi(i) .eq. 1)then
      do23022 j = 1,p
      ghib(j) = ghib(j) + a(j,i) * w(i)
23022 continue
23023 continue
      endif
      endif
      endif
23012 continue
23013 continue
      call dcopy(p,glob,1,aa(1,k+1),1)
      call dcopy(p,ghib,1,aa(1,k+2),1)
      yy(k+1) = big
      yy(k+2) = -big
      call dgemv('N',p,k+2,one-tau,aa,p,d,1,zero,rhs,1)
      call dscal(k+2,zero,wn,1)
      call daxpy(k+2,one-tau,u,1,wn,1)
      call rqfnb(k+2,p,aa,yy,rhs,d,u,beta,eps,wn,wp,nit(1,ir),info(ir))
      call dcopy(p,wp,1,b(1,ir),1)
      call dcopy(n,y,1,u,1)
      call dgemv('T',p,n,one,a,p,b(1,ir),1,one,u,1)
      sumbad = 0
      do23024 i = 1,n 
      if((u(i) .gt. 0) .and. slo(i) .eq. 1)then
      slo(i) = 0
      sumbad = sumbad + 1
      endif
      if((u(i) .lt. 0) .and. shi(i) .eq. 1)then
      shi(i) = 0
      sumbad = sumbad + 1
      endif
23024 continue
23025 continue
      if(sumbad .gt. 0)then
      if(sumbad .gt. 0.1 * nn)then
      nn = min(2 * nn, n)
      goto 23011
      endif
      else
      notopt = .false.
      endif
      goto 23010
      endif
23011 continue
      nit(4,ir) = ifix
      nit(5,ir) = ibad
      goto 23002
      endif
23003 continue
23000 continue
23001 continue
      return
      end
