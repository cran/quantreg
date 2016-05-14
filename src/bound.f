c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c  Function to obtain the step length
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine bound(x,dx,s,ds,z,dz,w,dw,n,beta,deltap,deltad)
c
      integer n
      double precision x(n),dx(n),s(n),ds(n),z(n),dz(n),w(n),dw(n)
      double precision deltap,deltad,dmin1,big,one,beta
      parameter (big = 1.0d20, one = 1.0d0)
      deltap = big
      deltad = big
      do i=1,n
         if(dx(i) .lt. 0) deltap = dmin1(deltap, -x(i)/dx(i))
         if(ds(i) .lt. 0) deltap = dmin1(deltap, -s(i)/ds(i))
         if(dz(i) .lt. 0) deltad = dmin1(deltad, -z(i)/dz(i))
         if(dw(i) .lt. 0) deltad = dmin1(deltad, -w(i)/dw(i))
      enddo
      deltap = dmin1(beta*deltap,one)
      deltad = dmin1(beta*deltad,one)
      return
      end
