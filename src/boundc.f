c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c  Function to obtain the step length
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine boundc(x1,dx1,x2,dx2,s,ds,z1,dz1,z2,dz2,w,dw,n1,n2,
     &                 beta,deltap,deltad)
c
      integer n1,n2
      double precision x1(n1),dx1(n1),x2(n2),dx2(n2),s(n1),ds(n1),
     &                 z1(n1),dz1(n1),z2(n2),dz2(n2),w(n1),dw(n1)
      double precision deltap,deltad,dmin1,big,one,beta
      parameter (big = 1.0d20, one = 1.0d0)
      deltap = big
      deltad = big
      do i = 1,n1
         if(dx1(i) .lt. 0) deltap = dmin1(deltap, -x1(i)/dx1(i))
         if(ds(i) .lt. 0) deltap = dmin1(deltap, -s(i)/ds(i))
         if(dz1(i) .lt. 0) deltad = dmin1(deltad, -z1(i)/dz1(i))
         if(dw(i) .lt. 0) deltad = dmin1(deltad, -w(i)/dw(i))
      enddo
      do i = 1,n2
         if(dx2(i) .lt. 0) deltap = dmin1(deltap, -x2(i)/dx2(i))
         if(dz2(i) .lt. 0) deltad = dmin1(deltad, -z2(i)/dz2(i))
      enddo
      deltap = dmin1(beta*deltap,one)
      deltad = dmin1(beta*deltad,one)
      return
      end
