      subroutine rqfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1,wn2,
&     wp,nit,info)
      integer n1,n2,p,info,nit(3)
      double precision a1(p,n1),a2(p,n2),y(n1),r(n2),rhs(p),d1(n1),d2(
&     n2),u(n1)
      double precision wn1(n1,9),wn2(n2,6),wp(p,p+3)
      double precision one,beta,eps
      data one/1.0d0/
      call lpfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1(1,1),wn2(1,
&     1),wn1(1,2),wp(1,1),wn1(1,3),wn2(1,2),wn1(1,4),wn1(1,5),wn2(1,3),
&     wn1(1,6), wp(1,2),wn1(1,7),wn2(1,4),wn1(1,8),wn1(1,9),wn2(1,5),
&     wn2(1,6),wp(1,3),wp(1,4),nit,info)
      return
      end
      subroutine lpfnc(n1,n2,p,a1,c1,a2,c2,b,d1,d2,u,beta,eps,x1,x2,s,y,
&     z1,z2,w,dx1,dx2,ds,dy,dz1,dz2,dw,dr1,dr2,r2, rhs,ada,nit,info)
      integer n1,p,i,info,nit(3),maxit
      double precision a1(p,n1),a2(p,n2),c1(n1),c2(n2),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz1,dxdz2,
&     dsdw
      double precision deltap,deltad,beta,eps,mu,gap,g
      double precision x1(n1),x2(n2),u(n1),s(n1),y(p),z1(n1),z2(n2),w(
&     n1)
      double precision d1(n1),d2(n2),rhs(p),ada(p,p)
      double precision dx1(n1),dx2(n2),ds(n1),dy(p),dz1(n1),dz2(n2),dw(
&     n1)
      double precision dr1(n1),dr2(n2),r2(n2)
      data zero /0.0d0/
      data one /1.0d0/
      data mone /-1.0d0/
      data big /1.0d+20/
      data maxit /50/
      nit(1)=0
      nit(2)=0
      nit(3)=n1
      call dgemv('N',p,n1,one,a1,p,c1,1,zero,y,1)
      do 23000 i=1,n1
      d1(i)=one
23000 continue
      do 23002 i=1,n2
      d2(i)=zero
      z2(i)=one
23002 continue
      call stepy2(n1,n2,p,a1,d1,a2,d2,y,ada,info)
      if(.not.(info .ne. 0))goto 23004
      return
23004 continue
      call dcopy(n1,c1,1,s,1)
      call dgemv('T',p,n1,mone,a1,p,y,1,one,s,1)
      do 23006 i=1,n1
      z1(i)=dmax1(s(i),zero)
      w(i)=dmax1(-s(i),zero)
      s(i)=u(i)-x1(i)
23006 continue
      gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
23008 if(.not.(gap .gt. eps .and. nit(1).lt.maxit))goto 23009
      nit(1)=nit(1)+1
      call dcopy(n2,c2,1,r2,1)
      call dgemv('T',p,n2,mone,a2,p,y,1,one,r2,1)
      call dcopy(p,b,1,dy,1)
      call dgemv('N',p,n1,mone,a1,p,x1,1,one,dy,1)
      call dgemv('N',p,n2,mone,a2,p,x2,1,one,dy,1)
      do 23010 i = 1,n1
      d1(i)=one/(z1(i)/x1(i) + w(i)/s(i))
      ds(i)=z1(i)-w(i)
      dz1(i)=d1(i)*ds(i)
23010 continue
      do 23012 i = 1,n2
      d2(i)=x2(i)/z2(i)
      dz2(i)=d2(i)*r2(i)
23012 continue
      call dgemv('N',p,n1,one,a1,p,dz1,1,one,dy,1)
      call dgemv('N',p,n2,one,a2,p,dz2,1,one,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy2(n1,n2,p,a1,d1,a2,d2,dy,ada,info)
      if(.not.(info .ne. 0))goto 23014
      return
23014 continue
      call dgemv('T',p,n1,one,a1,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do 23016 i=1,n1
      dx1(i)=d1(i)*ds(i)
      ds(i)=-dx1(i)
      dz1(i)=-z1(i)*(dx1(i)/x1(i) + one)
      dw(i)=-w(i)*(ds(i)/s(i) + one)
      if(.not.(dx1(i).lt.0))goto 23018
      deltap=dmin1(deltap,-x1(i)/dx1(i))
23018 continue
      if(.not.(ds(i).lt.0))goto 23020
      deltap=dmin1(deltap,-s(i)/ds(i))
23020 continue
      if(.not.(dz1(i).lt.0))goto 23022
      deltad=dmin1(deltad,-z1(i)/dz1(i))
23022 continue
      if(.not.(dw(i).lt.0))goto 23024
      deltad=dmin1(deltad,-w(i)/dw(i))
23024 continue
23016 continue
      call dcopy(n2,r2,1,dx2,1)
      call dgemv('T',p,n2,one,a2,p,dy,1,mone,dx2,1)
      do 23026 i=1,n2
      dx2(i)=d2(i)*dx2(i)
      dz2(i)=-z2(i)*(dx2(i)/x2(i) + one)
      if(.not.(dx2(i).lt.0))goto 23028
      deltap=dmin1(deltap,-x2(i)/dx2(i))
23028 continue
      if(.not.(dz2(i).lt.0))goto 23030
      deltad=dmin1(deltad,-z2(i)/dz2(i))
23030 continue
23026 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(.not.(min(deltap,deltad) .lt. one))goto 23032
      nit(2)=nit(2)+1
      mu = ddot(n1,x1,1,z1,1)+ddot(n2,x2,1,z2,1)+ddot(n1,s,1,w,1)
      g = mu + deltap*ddot(n1,dx1,1,z1,1)+deltad*ddot(n1,dz1,1,x1,1)+
&     deltap*deltad*ddot(n1,dz1,1,dx1,1)+deltap*ddot(n2,dx2,1,z2,1)+
&     deltad*ddot(n2,dz2,1,x2,1)+deltap*deltad*ddot(n2,dz2,1,dx2,1)+
&     deltap*ddot(n1,ds,1,w,1)+deltad*ddot(n1,dw,1,s,1) +deltap*deltad*
&     ddot(n1,ds,1,dw,1)
      mu = mu * ((g/mu)**3) /(dfloat(2*n1)+dfloat(n2))
      do 23034 i=1,n1
      dsdw = ds(i)*dw(i)
      dr1(i)=d1(i)*(mu*(one/s(i)-one/x1(i))+dx1(i)*dz1(i)/x1(i)-dsdw/s(
&     i))
23034 continue
      do 23036 i=1,n2
      dr2(i)=d2(i)*(dx2(i)*dz2(i)/x2(i)-mu/x2(i))
23036 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n1,one,a1,p,dr1,1,one,dy,1)
      call dgemv('N',p,n2,one,a2,p,dr2,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call dgemv('T',p,n1,one,a1,p,dy,1,zero,u,1)
      deltap=big
      deltad=big
      do 23038 i=1,n1
      dsdw = ds(i)*dw(i)
      dxdz1 = dx1(i)*dz1(i)
      dx1(i) = d1(i)*(u(i)-z1(i)+w(i))-dr1(i)
      ds(i) = -dx1(i)
      dz1(i) = -z1(i)+(mu - z1(i)*dx1(i) - dxdz1)/x1(i)
      dw(i) = -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
      if(.not.(dx1(i).lt.0))goto 23040
      deltap=dmin1(deltap,-x1(i)/dx1(i))
23040 continue
      if(.not.(ds(i).lt.0))goto 23042
      deltap=dmin1(deltap,-s(i)/ds(i))
23042 continue
      if(.not.(dz1(i).lt.0))goto 23044
      deltad=dmin1(deltad,-z1(i)/dz1(i))
23044 continue
      if(.not.(dw(i).lt.0))goto 23046
      deltad=dmin1(deltad,-w(i)/dw(i))
23046 continue
23038 continue
      call dgemv('T',p,n2,one,a2,p,dy,1,zero,u,1)
      do 23048 i=1,n2
      dxdz2 = dx2(i)*dz2(i)
      dx2(i) = d2(i)*(u(i)-r2(i))-dr2(i)
      dz2(i) = -z2(i)+(mu - z2(i)*dx2(i) - dxdz2)/x2(i)
      if(.not.(dx2(i).lt.0))goto 23050
      deltap=dmin1(deltap,-x2(i)/dx2(i))
23050 continue
      if(.not.(dz2(i).lt.0))goto 23052
      deltad=dmin1(deltad,-z2(i)/dz2(i))
23052 continue
23048 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
23032 continue
      call daxpy(n1,deltap,dx1,1,x1,1)
      call daxpy(n2,deltap,dx2,1,x2,1)
      call daxpy(n1,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n1,deltad,dz1,1,z1,1)
      call daxpy(n2,deltad,dz2,1,z2,1)
      call daxpy(n1,deltad,dw,1,w,1)
      gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
      goto 23008
23009 continue
      call daxpy(n1,mone,w,1,z1,1)
      call dswap(n1,z1,1,x1,1)
      return
      end
      subroutine stepy2(n1,n2,p,a1,d1,a2,d2,b,ada,info)
      integer n1,n2,p,i,j,k,info
      double precision a1(p,n1),a2(p,n2),b(p),d1(n1),d2(n2),ada(p,p),
&     zero
      data zero/0.0d0/
      do 23054 j=1,p
      do 23056 k=1,p
      ada(j,k)=zero
23056 continue
23054 continue
      do 23058 i=1,n1
      call dsyr('U',p,d1(i),a1(1,i),1,ada,p)
23058 continue
      do 23060 i=1,n2
      call dsyr('U',p,d2(i),a2(1,i),1,ada,p)
23060 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end
