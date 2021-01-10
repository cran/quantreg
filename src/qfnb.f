C Output from Public domain Ratfor, version 1.05
      subroutine qfnb(n,p,m,a,y,t,r,d,u,wn,wp,b,nit, info)
      integer n,p,m,nit(3),info
      double precision a(p,n), y(n), t(m), b(p,m), r(p)
      double precision d(n), u(n), wn(n,9), wp(p,p+3)
      double precision zero, one, eps, beta
      parameter( zero = 0.0d0)
      parameter( one = 1.0d0)
      parameter( beta = 0.99995d0)
      parameter( eps = 1.0d-6)
      do23000 i = 1,m
      call dgemv('N',p,n,one-t(i),a,p,d,1,zero,r,1)
      call dscal(n,zero,wn,1)
      call daxpy(n,one-t(i),u,1,wn,1)
      call rqfnb(n,p,a,y,r,d,u,beta,eps,wn,wp,nit,info)
      if(info .ne. 0)then
      goto 23001
      endif
      do23004 j = 1,n
      u(j) = one
      d(j) = one
23004 continue
23005 continue
      call dcopy(p,wp,1,b(1,i),1)
23000 continue
23001 continue
      return
      end
