      subroutine rqfnb(n,p,a,y,rhs,d,u,beta,eps,wn,wp,nit,info)
      integer n,p,info,nit(3)
      double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,9),wp(p,p+3)
      double precision one,beta,eps
      data one/1.0d0/
      call lpfnb(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2), wp(1,1),wn(1,
&     3),wn(1,4),wn(1,5),wn(1,6), wp(1,2),wn(1,7),wn(1,8),wn(1,9),wp(1,
&     3),wp(1,4),nit,info)
      return
      end
      subroutine lpfnb(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w,dx,ds,dy,dz,dw,
&     dr,rhs,ada,nit,info)
      integer n,p,pp,i,info,nit(3),maxit
      double precision a(p,n),c(n),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz,dsdw
      double precision deltap,deltad,beta,eps,mu,gap,g
      double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,
&     p)
      double precision dx(n),ds(n),dy(p),dz(n),dw(n),dr(n)
      data zero /0.0d0/
      data one /1.0d0/
      data mone /-1.0d0/
      data big /1.0d+20/
      data maxit /50/
      nit(1)=0
      nit(2)=0
      nit(3)=n
      pp=p*p
      call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
      do 23000 i=1,n
      d(i)=one
23000 continue
      call stepy(n,p,a,d,y,ada,info)
      if(.not.(info .ne. 0))goto 23002
      return
23002 continue
      call dcopy(n,c,1,s,1)
      call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
      do 23004 i=1,n
      z(i)=dmax1(s(i),zero)
      w(i)=dmax1(-s(i),zero)
      s(i)=u(i)-x(i)
23004 continue
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
23006 if(.not.(gap .gt. eps .and. nit(1).lt.maxit))goto 23007
      nit(1)=nit(1)+1
      do 23008 i = 1,n
      d(i) = one/(z(i)/x(i) + w(i)/s(i))
      ds(i)=z(i)-w(i)
      dz(i)=d(i)*ds(i)
23008 continue
      call dgemv('N',p,n,one,a,p,dz,1,zero,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy(n,p,a,d,dy,ada,info)
      if(.not.(info .ne. 0))goto 23010
      return
23010 continue
      call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do 23012 i=1,n
      dx(i)=d(i)*ds(i)
      ds(i)=-dx(i)
      dz(i)=-z(i)*(dx(i)/x(i) + one)
      dw(i)=-w(i)*(ds(i)/s(i) + one)
      if(.not.(dx(i).lt.0))goto 23014
      deltap=dmin1(deltap,-x(i)/dx(i))
23014 continue
      if(.not.(ds(i).lt.0))goto 23016
      deltap=dmin1(deltap,-s(i)/ds(i))
23016 continue
      if(.not.(dz(i).lt.0))goto 23018
      deltad=dmin1(deltad,-z(i)/dz(i))
23018 continue
      if(.not.(dw(i).lt.0))goto 23020
      deltad=dmin1(deltad,-w(i)/dw(i))
23020 continue
23012 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(.not.(min(deltap,deltad) .lt. one))goto 23022
      nit(2)=nit(2)+1
      mu = ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
      g = mu + deltap*ddot(n,dx,1,z,1)+deltad*ddot(n,dz,1,x,1) +deltap*
&     deltad*ddot(n,dz,1,dx,1)+deltap*ddot(n,ds,1,w,1)+deltad*ddot(n,dw,
&     1,s,1) +deltap*deltad*ddot(n,ds,1,dw,1)
      mu = mu * ((g/mu)**3) /dfloat(2*n)
      do 23024 i=1,n
      dr(i)=d(i)*(mu*(1/s(i)-1/x(i))+dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
23024 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n,one,a,p,dr,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call dgemv('T',p,n,one,a,p,dy,1,zero,u,1)
      deltap=big
      deltad=big
      do 23026 i=1,n
      dxdz = dx(i)*dz(i)
      dsdw = ds(i)*dw(i)
      dx(i)= d(i)*(u(i)-z(i)+w(i))-dr(i)
      ds(i)= -dx(i)
      dz(i)= -z(i)+(mu - z(i)*dx(i) - dxdz)/x(i)
      dw(i)= -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
      if(.not.(dx(i).lt.0))goto 23028
      deltap=dmin1(deltap,-x(i)/dx(i))
23028 continue
      if(.not.(ds(i).lt.0))goto 23030
      deltap=dmin1(deltap,-s(i)/ds(i))
23030 continue
      if(.not.(dz(i).lt.0))goto 23032
      deltad=dmin1(deltad,-z(i)/dz(i))
23032 continue
      if(.not.(dw(i).lt.0))goto 23034
      deltad=dmin1(deltad,-w(i)/dw(i))
23034 continue
23026 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
23022 continue
      call daxpy(n,deltap,dx,1,x,1)
      call daxpy(n,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n,deltad,dz,1,z,1)
      call daxpy(n,deltad,dw,1,w,1)
      gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
      goto 23006
23007 continue
      call daxpy(n,mone,w,1,z,1)
      call dswap(n,z,1,x,1)
      return
      end
      subroutine stepy(n,p,a,d,b,ada,info)
      integer n,p,pp,i,info
      double precision a(p,n),b(p),d(n),ada(p,p),zero
      data zero/0.0d0/
      pp=p*p
      do 23036 j=1,p
      do 23038 k=1,p
      ada(j,k)=zero
23038 continue
23036 continue
      do 23040 i=1,n
      call dsyr('U',p,d(i),a(1,i),1,ada,p)
23040 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end
