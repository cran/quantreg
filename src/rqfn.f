C Output from Public domain Ratfor, version 1.0
      subroutine rqfn(n,p,a,y,rhs,d,u,beta,eps,wn,wp,aa,nit,info)
      integer n,p,info,nit(3)
      double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,10),wp(p,p+3),a
     *a(p,p)
      double precision one,beta,eps
      parameter( one = 1.0d0)
      call fna(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2), wp(1,1),wn(1,3)
     *,wn(1,4),wn(1,5), wn(1,6), wp(1,2),wn(1,7),wn(1,8),wn(1,9),wn(1,10
     *),wp(1,3), wp(1,4),aa,nit,info)
      return
      end
      subroutine fna(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w, dx,ds,dy,dz,dw,ds
     *dw,dxdz,rhs,ada,aa,nit,info)
      integer n,p,pp,i,info,nit(3)
      double precision a(p,n),c(n),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dasum
      double precision deltap,deltad,beta,eps,cx,by,uw,uz,mu,mua,acomp,r
     *dg,g
      double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p
     *)
      double precision aa(p,p),dx(n),ds(n),dy(p),dz(n),dw(n),dxdz(n),dsd
     *w(n)
      parameter( zero = 0.0d0)
      parameter( half = 0.5d0)
      parameter( one = 1.0d0)
      parameter( mone = -1.0d0)
      parameter( big = 1.0d+20)
      nit(1)=0
      nit(2)=0
      nit(3)=n
      pp=p*p
      call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
      call stepy(n,p,a,d,y,aa,info)
      if(info .ne. 0)then
      return
      endif
      do23002 i=1,p
      do23004 j=1,p
      ada(i,j)=zero
23004 continue
23005 continue
      ada(i,i)=one
23002 continue
23003 continue
      call dtrtrs('U','T','N',p,p,aa,p,ada,p,info)
      call dcopy(pp,ada,1,aa,1)
      call dcopy(n,c,1,s,1)
      call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
      do23006 i=1,n
      d(i)=one
      if(dabs(s(i)) .lt. eps)then
      z(i) = dmax1( s(i),zero) + eps
      w(i) = dmax1(-s(i),zero) + eps
      else
      z(i) = dmax1( s(i),zero)
      w(i) = dmax1(-s(i),zero)
      endif
      s(i)=u(i)-x(i)
23006 continue
23007 continue
      cx = ddot(n,c,1,x,1)
      by = ddot(p,b,1,y,1)
      uw = dasum(n,w,1)
      uz = dasum(n,z,1)
      rdg = (cx - by + uw)
23010 if(rdg .gt. eps)then
      nit(1)=nit(1)+1
      do23012 i =1,n
      d(i) = one/(z(i)/x(i) + w(i)/s(i))
      ds(i)=z(i)-w(i)
      dx(i)=d(i)*ds(i)
23012 continue
23013 continue
      call dgemv('N',p,n,one,a,p,dx,1,zero,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy(n,p,a,d,dy,ada,info)
      if(info .ne. 0)then
      return
      endif
      call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do23016 i=1,n
      dx(i)=d(i)*ds(i)
      ds(i)=-dx(i)
      dz(i)=-z(i)*(dx(i)/x(i) + one)
      dw(i)=w(i)*(dx(i)/s(i) - one)
      dxdz(i)=dx(i)*dz(i)
      dsdw(i)=ds(i)*dw(i)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23016 continue
23017 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(deltap*deltad.lt.one)then
      nit(2)=nit(2)+1
      acomp=ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
      g=acomp+deltap*ddot(n,dx,1,z,1)+ deltad*ddot(n,dz,1,x,1)+ deltap*d
     *eltad*ddot(n,dz,1,dx,1)+ deltap*ddot(n,ds,1,w,1)+ deltad*ddot(n,dw
     *,1,s,1)+ deltap*deltad*ddot(n,ds,1,dw,1)
      mu=acomp/dfloat(2*n)
      mua=g/dfloat(2*n)
      mu=mu*(mua/mu)**3
      do23028 i=1,n
      dz(i)=d(i)*(mu*(1/s(i)-1/x(i))+ dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
23028 continue
23029 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call daxpy(p,mone,dy,1,rhs,1)
      call dgemv('T',p,n,one,a,p,rhs,1,zero,dw,1)
      deltap=big
      deltad=big
      do23030 i=1,n
      dx(i)=dx(i)-dz(i)-d(i)*dw(i)
      ds(i)=-dx(i)
      dz(i)=mu/x(i) - z(i)*dx(i)/x(i) - z(i) - dxdz(i)/x(i)
      dw(i)=mu/s(i) - w(i)*ds(i)/s(i) - w(i) - dsdw(i)/s(i)
      if(dx(i).lt.0)then
      deltap=dmin1(deltap,-x(i)/dx(i))
      else
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz(i).lt.0)then
      deltad=dmin1(deltad,-z(i)/dz(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23030 continue
23031 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      endif
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
      goto 23010
      endif
23011 continue
      call daxpy(n,mone,w,1,z,1)
      call dswap(n,z,1,x,1)
      return
      end
