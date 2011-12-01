C Output from Public domain Ratfor, version 1.0
      subroutine brutpow(n,p,m,h,a,b,c,x,tau,u,xh,d,jminz,nflag)
      integer n,p,m
      double precision x(p),a(n,p),b(n),c(n)
      double precision u(p,p),d(p),xh(p)
      double precision zero, one,tau,pow,minz,z
      integer h(p,m),k,findk,jminz,nflag
      parameter(zero = 0.0d0, one = 1.d0)
      minz = pow(n,p,x,a,b,c,tau)
      do23000 j = 2,m 
      k = findk(p,h(1,j),h(1,j-1))
      if(k .eq. 0)then
      nflag = 4
      return
      endif
      call pivot(n,p,h(1,j-1),h(k,j),h(k,j-1),a,u,d,xh,nflag)
      if(nflag .gt. 0)then
      return
      endif
      do23006 i = 1,p
      xh(i) = b(h(i,j))
23006 continue
23007 continue
      call dgemv('N',p,p,one,u,p,xh,1,zero,x,1)
      z = pow(n,p,x,a,b,c,tau)
      if(z .lt. minz)then
      minz = z
      jminz = j
      endif
23000 continue
23001 continue
      return
      end
      integer function findk(p,h,g)
      integer p,k,h(p),g(p)
      findk = 0
      do23010 k = 1,p
      if(h(k) .ne. g(k))then
      findk = k
      goto 23011
      endif
23010 continue
23011 continue
      return
      end
