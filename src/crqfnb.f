C Output from Public domain Ratfor, version 1.0
      subroutine crqfnb(n,p,a1,c1,n1,x,y,c,b,g,m,r,s,d,u,wn,wp,info)
      integer n,p,n1,m,info,nit(3)
      double precision a1(p,n1),c1(n),x(n,p),y(n),c(n),b(p,m),g(m)
      double precision wn(n,9),wp(p,p+3),r(p),s(n),d(n),u(n)
      double precision zero,half,one,beta,eps,dh
      parameter( zero = 0.0d0)
      parameter( half = 0.5d0)
      parameter( one = 1.0d0)
      parameter( beta = 0.99995d0)
      parameter( eps = 1.0d-8)
      do23000 k = 2,m 
      dh = -log(one - g(k)) + log(one - g(k-1))
      do23002 i = 1,n 
      u(i) = one
      wn(i,1) = half
      if(d(i) .ge. zero)then
      s(i) = s(i) + dh
      endif
      d(i) = c(i) - s(i)
23002 continue
23003 continue
      call dgemv('T',n,p,one,x,n,d,1,zero,r,1)
      call rqfnb(n1,p,a1,c1,r,d,u,beta,eps,wn,wp,nit,info)
      if(info .ne. 0)then
      goto 23001
      endif
      call dcopy(p,wp,1,b(1,k-1),1)
      call dcopy(n,y,1,d,1)
      call dgemv('N',n,p,one,x,n,b(1,k-1),1,one,d,1)
23000 continue
23001 continue
      m = k-1
      return
      end
