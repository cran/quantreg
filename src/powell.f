C Output from Public domain Ratfor, version 1.0
      subroutine powell(n,p,p2,a,b,c,x,tau,h,f,u,s,g,d,xh,maxit,nflag)
      integer n,p,p2
      double precision x(p),a(n,p),b(n),c(n)
      double precision f(n),u(p,p),s(n),g(p2),d(p),xh(p)
      double precision zero, one, mone, step, tau,pow
      integer h(p),hin,hout,k,it,inset,maxit,nflag
      parameter(zero = 0.0d0, one = 1.d0, mone = -1.d0)
      it = 0
23000 continue
      it = it + 1
      if(it .gt. 1)then
      call pivot(n,p,h,hin,hout,a,u,d,xh,nflag)
      endif
      if(nflag .gt. 0)then
      nflag = nflag + 2
      return
      endif
      do23007 i = 1,p
      xh(i) = b(h(i))
23007 continue
23008 continue
      call dgemv('N',p,p,one,u,p,xh,1,zero,x,1)
      call dgemv('N',n,p,one,a,n,x,1,zero,f,1)
      do23009 i = 1,n
      if(inset(p,i,h) .gt. 0 .or. f(i) .gt. c(i))then
      s(i) = zero
      else
      if(b(i) .lt. f(i))then
      s(i) = one - tau
      else
      s(i) = - tau
      endif
      endif
23009 continue
23010 continue
      call dgemv('T',n,p,one,a,n,s,1,zero,xh,1)
      call dgemv('T',p,p,one,u,p,xh,1,zero,g,1)
      do23015 i = 1,p 
      if(f(h(i)) .lt. c(h(i)))then
      if(b(h(i)) .lt. c(h(i)))then
      g(i + p) = - g(i) + one - tau
      else
      g(i + p) = - g(i) - tau
      endif
      else
      g(i + p) = - g(i) + tau
      endif
      g(i) = g(i) + one - tau
23015 continue
23016 continue
      k = idmin(p2,g,1)
      if(g(k) .ge. 0 .or. it .gt. maxit)then
      goto 23002
      endif
      call dscal(p,zero,d,1)
      if(k .le. p)then
      call daxpy(p,one,u(1,k),1,d,1)
      else
      k = k - p
      call daxpy(p,mone,u(1,k),1,d,1)
      endif
      call dgemv('N',n,p,one,a,n,d,1,zero,s,1)
      do23025 i = 1,n 
      call dcopy(p,x,1,xh,1)
      step = (b(i) - f(i))/s(i)
      call daxpy(p,step,d,1,xh,1)
      s(i) = pow(n,p,xh,a,b,c,tau)
23025 continue
23026 continue
      hin = idmin(n,s,1)
      if(inset(p,hin,h) .gt. 0)then
      nflag = 2
      goto 23002
      endif
      hout = h(k)
23001 goto 23000
23002 continue
      if(it .gt. maxit)then
      nflag = 1
      endif
      return
      end
      subroutine pivot(n,p,h,hin,hout,a,b,u,v,eflag)
      integer n,p,h(p),hin,hout,inset,k,eflag
      double precision a(n,p),b(p,p),u(p),v(p)
      double precision zero,one
      parameter(zero = 0.d0, one = 1.d0)
      eflag = 0
      k = inset(p,hout,h)
      if(k .eq. 0)then
      eflag = 1
      return
      endif
      if(inset(p,hin,h) .gt. 0)then
      eflag = 2
      return
      endif
      if(hin .lt. 1 .or. hin .gt. n)then
      eflag = 3
      return
      endif
      call dcopy(p,a(hin,1),n,v,1)
      call dgemv('T',p,p,one,b,p,v,1,zero,u,1)
      call dcopy(p,b(1,k),1,v,1)
      do23037 j = 1,p
      do23039 i = 1,p
      if(j .eq. k)then
      b(i,j) = b(i,j)/u(k)
      else
      b(i,j) = b(i,j) - (u(j)/u(k)) * v(i)
      endif
23039 continue
23040 continue
23037 continue
23038 continue
      h(k) = hin
      return
      end
      integer function inset(p,k,h)
      integer p,k,h(p)
      do23043 inset = 1,p
      if(h(inset) .eq. k)then
      return
      endif
23043 continue
23044 continue
      inset = 0
      return
      end
      double precision function pow(n,p,x,a,b,c,tau)
      integer n,p
      double precision x(p),a(n,p),b(n),c(n)
      double precision tau,u,zero,rho,fit,ddot
      parameter(zero= 0.d0)
      pow = zero
      do23047 i = 1,n
      fit = ddot(p,a(i,1),n,x,1)
      u = b(i) - min(fit,c(i))
      pow = pow + rho(u, tau)
23047 continue
23048 continue
      return
      end
      double precision function rho(u,tau)
      double precision u,tau,one
      parameter(one = 1.d0)
      if(u .lt. 0)then
      rho = u * (tau - one)
      else
      rho = u * tau
      endif
      return
      end
