      subroutine rls(n,p,x,y,b,a,ax)
      integer n,p
      double precision x(p,n),y(n),b(p,n),a(p,p),ax(p)
      double precision zero,one,mone,f,r,ddot
      parameter( one = 1.d0)
      parameter( mone = -1.d0)
      parameter( zero = 0.d0)
      do 23000 i = (p+1),n 
      call dgemv('N',p,p,one,a,p,x(1,i),1,zero,ax,1)
      f = one + ddot(p,x(1,i),1,ax,1)
      r = (y(i)-ddot(p,x(1,i),1,b(1,i-1),1))/f
      call daxpy(p,one,b(1,i-1),1,b(1,i),1)
      call daxpy(p,r,ax,1,b(1,i),1)
      call dger(p,p,mone/f,ax,1,ax,1,a,p)
23000 continue
      return
      end
