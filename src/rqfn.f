      subroutine rqfn(n,p,a,y,rhs,d,beta,eps,tau,wn,wp,aa,nit,info)
      integer n,p,info,nit(3)
      double precision a(p,n),y(n),rhs(p),d(n),beta,eps,wn(n,10),wp(p,p+
     &3),aa(p,p)
      double precision mone,one,tau,ddot
      data one/1.0d0/
      data mone/-1.0d0/
      do 23000 i=1,n
      d(i)=one
      wn(i,1)=one-tau
23000 continue
      do 23002 i=1,p
      rhs(i)=(one-tau)*ddot(n,d,1,a(i,1),p)
23002 continue
      call dscal(n,mone,y,1)
      call fna(n,p,a,y,rhs,d,beta,eps,wn(1,1),wn(1,2),wp(1,1),wn(1,3),
     &wn(1,4),wn(1,5), wn(1,6),wp(1,2),wn(1,7),wn(1,8),wn(1,9),wn(1,10),
     &wp(1,3), wp(1,4),aa,nit,info)
      call dscal(p,mone,wp,1)
      call dscal(n,mone,wn,1)
      call dscal(n,mone,y,1)
      return
      end
      subroutine fna(n,p,a,c,b,d,beta,eps,x,s,y,z,w,dx,ds,dy,dz,dw,dsdw,
     &dxdz,rhs,ada,aa,nit,info)
      integer n,p,pp,i,info,nit(3)
      double precision a(p,n),c(n),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dasum
      double precision deltap,deltad,beta,eps,cx,by,uw,uz,mu,mua,acomp,
     &rdg,g
      double precision x(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p),aa(
     &p,p)
      double precision dx(n),ds(n),dy(p),dz(n),dw(n),dxdz(n),dsdw(n)
      data zero /0.0d0/
      data one /1.0d0/
      data mone /-1.0d0/
      data big /1.0d+20/
      nit(1)=0
      nit(2)=0
      nit(3)=n
      pp=p*p
      call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
      call stepy(n,p,a,d,y,aa,info)
      if(.not.(info .ne. 0))goto 23004
      return
23004 continue
      do 23006 i=1,p
      do 23008 j=1,p
      ada(i,j)=zero
23008 continue
      ada(i,i)=one
23006 continue
      call dtrtrs('U','T','N',p,p,aa,p,ada,p,info)
      call dcopy(pp,ada,1,aa,1)
      call dcopy(n,c,1,s,1)
      call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
      do 23010 i=1,n
      d(i)=one
      z(i)=dmax1(s(i),zero)
      w(i)=dmax1(-s(i),zero)
      s(i)=one-x(i)
23010 continue
      cx = ddot(n,c,1,x,1)
      by = ddot(p,b,1,y,1)
      uw = dasum(n,w,1)
      uz = dasum(n,z,1)
      rdg = (cx - by + uw)
23012 if(.not.(rdg .gt. eps))goto 23013
      nit(1)=nit(1)+1
      do 23014 i =1,n
      d(i) = one/(z(i)/x(i) + w(i)/s(i))
      ds(i)=z(i)-w(i)
      dx(i)=d(i)*ds(i)
23014 continue
      call dgemv('N',p,n,one,a,p,dx,1,zero,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy(n,p,a,d,dy,ada,info)
      if(.not.(info .ne. 0))goto 23016
      return
23016 continue
      call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do 23018 i=1,n
      dx(i)=d(i)*ds(i)
      ds(i)=-dx(i)
      dz(i)=-z(i)*(dx(i)/x(i) + one)
      dw(i)=w(i)*(dx(i)/s(i) - one)
      dxdz(i)=dx(i)*dz(i)
      dsdw(i)=ds(i)*dw(i)
      if(.not.(dx(i).lt.0))goto 23020
      deltap=dmin1(deltap,-x(i)/dx(i))
23020 continue
      if(.not.(ds(i).lt.0))goto 23022
      deltap=dmin1(deltap,-s(i)/ds(i))
23022 continue
      if(.not.(dz(i).lt.0))goto 23024
      deltad=dmin1(deltad,-z(i)/dz(i))
23024 continue
      if(.not.(dw(i).lt.0))goto 23026
      deltad=dmin1(deltad,-w(i)/dw(i))
23026 continue
23018 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(.not.(deltap*deltad.lt.one))goto 23028
      nit(2)=nit(2)+1
      acomp=ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
      g=acomp+deltap*ddot(n,dx,1,z,1)+deltad*ddot(n,dz,1,x,1)+deltap*
     &deltad*ddot(n,dz,1,dx,1)+deltap*ddot(n,ds,1,w,1)+deltad*ddot(n,dw,
     &1,s,1)+deltap*deltad*ddot(n,ds,1,dw,1)
      mu=acomp/dfloat(2*n)
      mua=g/dfloat(2*n)
      mu=mu*(mua/mu)**3
      do 23030 i=1,n
      dz(i)=d(i)*(mu*(1/s(i)-1/x(i))+dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
23030 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call daxpy(p,mone,dy,1,rhs,1)
      call dgemv('T',p,n,one,a,p,rhs,1,zero,dw,1)
      deltap=big
      deltad=big
      do 23032 i=1,n
      dx(i)=dx(i)-dz(i)-d(i)*dw(i)
      ds(i)=-dx(i)
      dz(i)=mu/x(i) - z(i)*dx(i)/x(i) - z(i) - dxdz(i)/x(i)
      dw(i)=mu/s(i) - w(i)*ds(i)/s(i) - w(i) - dsdw(i)/s(i)
      if(.not.(dx(i).lt.0))goto 23034
      deltap=dmin1(deltap,-x(i)/dx(i))
      goto 23035
23034 continue
      deltap=dmin1(deltap,-s(i)/ds(i))
23035 continue
      if(.not.(dz(i).lt.0))goto 23036
      deltad=dmin1(deltad,-z(i)/dz(i))
23036 continue
      if(.not.(dw(i).lt.0))goto 23038
      deltad=dmin1(deltad,-w(i)/dw(i))
23038 continue
23032 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
23028 continue
      call daxpy(n,deltap,dx,1,x,1)
      call daxpy(n,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n,deltad,dz,1,z,1)
      call daxpy(n,deltad,dw,1,w,1)
      cx=ddot(n,c,1,x,1)
      by=ddot(p,b,1,y,1)
      uw = dasum(n,w,1)
      uz = dasum(n,z,1)
      rdg=(cx-by+uw)
      goto 23012
23013 continue
      call daxpy(n,mone,w,1,z,1)
      call dswap(n,z,1,x,1)
      return
      end
      subroutine stepy(n,p,a,d,b,ada,info)
      integer n,p,pp,i,info
      double precision a(p,n),b(p),d(n),ada(p,p),zero
      data zero/0.0d0/
      pp=p*p
      do 23040 j=1,p
      do 23042 k=1,p
      ada(j,k)=zero
23042 continue
23040 continue
      do 23044 i=1,n
      call dsyr('U',p,d(i),a(1,i),1,ada,p)
23044 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end
