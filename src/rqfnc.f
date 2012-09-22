C Output from Public domain Ratfor, version 1.0
      subroutine rqfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1,wn2,wp
     *,nit,info)
      integer n1,n2,p,info,nit(3)
      double precision a1(p,n1),a2(p,n2),y(n1),r(n2),rhs(p),d1(n1),d2(n2
     *),u(n1)
      double precision wn1(n1,9),wn2(n2,6),wp(p,p+3)
      double precision one,beta,eps
      parameter(one = 1.0d0)
      call lpfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1(1,1),wn2(1,1
     *),wn1(1,2), wp(1,1),wn1(1,3),wn2(1,2),wn1(1,4),wn1(1,5),wn2(1,3),w
     *n1(1,6), wp(1,2),wn1(1,7),wn2(1,4),wn1(1,8),wn1(1,9),wn2(1,5),wn2(
     *1,6), wp(1,3),wp(1,4),nit,info)
      return
      end
      subroutine lpfnc(n1,n2,p,a1,c1,a2,c2,b,d1,d2,u,beta,eps,x1,x2,s, y
     *,z1,z2,w,dx1,dx2,ds,dy,dz1,dz2,dw,dr1,dr2,r2, rhs,ada,nit,info)
      integer n1,p,i,info,nit(3),maxit
      double precision a1(p,n1),a2(p,n2),c1(n1),c2(n2),b(p)
      double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz1,dxdz2,ds
     *dw
      double precision deltap,deltad,beta,eps,mu,gap,g
      double precision x1(n1),x2(n2),u(n1),s(n1),y(p),z1(n1),z2(n2),w(n1
     *)
      double precision d1(n1),d2(n2),rhs(p),ada(p,p)
      double precision dx1(n1),dx2(n2),ds(n1),dy(p),dz1(n1),dz2(n2),dw(n
     *1)
      double precision dr1(n1),dr2(n2),r2(n2)
      parameter(zero = 0.0d0)
      parameter(one = 1.0d0)
      parameter(mone = -1.0d0)
      parameter(big = 1.0d+20)
      parameter(maxit = 500)
      nit(1)=0
      nit(2)=0
      nit(3)=n1
      call dgemv('N',p,n1,one,a1,p,c1,1,zero,y,1)
      do23000 i=1,n1
      d1(i)=one
23000 continue
23001 continue
      do23002 i=1,n2
      d2(i)=zero
      z2(i)=one
23002 continue
23003 continue
      call stepy2(n1,n2,p,a1,d1,a2,d2,y,ada,info)
      if(info .ne. 0)then
      return
      endif
      call dcopy(n1,c1,1,s,1)
      call dgemv('T',p,n1,mone,a1,p,y,1,one,s,1)
      do23006 i=1,n1
      if(dabs(s(i)) .lt. eps)then
      z1(i)=dmax1(s(i),zero)+eps
      w(i)=dmax1(-s(i),zero)+eps
      else
      z1(i)=dmax1(s(i),zero)
      w(i)=dmax1(-s(i),zero)
      endif
      s(i)=u(i)-x1(i)
23006 continue
23007 continue
      gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
23010 if(gap .gt. eps .and. nit(1).lt.maxit)then
      nit(1)=nit(1)+1
      call dcopy(n2,c2,1,r2,1)
      call dgemv('T',p,n2,mone,a2,p,y,1,one,r2,1)
      call dcopy(p,b,1,dy,1)
      call dgemv('N',p,n1,mone,a1,p,x1,1,one,dy,1)
      call dgemv('N',p,n2,mone,a2,p,x2,1,one,dy,1)
      do23012 i = 1,n1
      d1(i)=one/(z1(i)/x1(i) + w(i)/s(i))
      ds(i)=z1(i)-w(i)
      dz1(i)=d1(i)*ds(i)
23012 continue
23013 continue
      do23014 i = 1,n2
      d2(i)=x2(i)/z2(i)
      dz2(i)=d2(i)*r2(i)
23014 continue
23015 continue
      call dgemv('N',p,n1,one,a1,p,dz1,1,one,dy,1)
      call dgemv('N',p,n2,one,a2,p,dz2,1,one,dy,1)
      call dcopy(p,dy,1,rhs,1)
      call stepy2(n1,n2,p,a1,d1,a2,d2,dy,ada,info)
      if(info .ne. 0)then
      return
      endif
      call dgemv('T',p,n1,one,a1,p,dy,1,mone,ds,1)
      deltap=big
      deltad=big
      do23018 i=1,n1
      dx1(i)=d1(i)*ds(i)
      ds(i)=-dx1(i)
      dz1(i)=-z1(i)*(dx1(i)/x1(i) + one)
      dw(i)=-w(i)*(ds(i)/s(i) + one)
      if(dx1(i).lt.0)then
      deltap=dmin1(deltap,-x1(i)/dx1(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz1(i).lt.0)then
      deltad=dmin1(deltad,-z1(i)/dz1(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23018 continue
23019 continue
      call dcopy(n2,r2,1,dx2,1)
      call dgemv('T',p,n2,one,a2,p,dy,1,mone,dx2,1)
      do23028 i=1,n2
      dx2(i)=d2(i)*dx2(i)
      dz2(i)=-z2(i)*(dx2(i)/x2(i) + one)
      if(dx2(i).lt.0)then
      deltap=dmin1(deltap,-x2(i)/dx2(i))
      endif
      if(dz2(i).lt.0)then
      deltad=dmin1(deltad,-z2(i)/dz2(i))
      endif
23028 continue
23029 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      if(min(deltap,deltad) .lt. one)then
      nit(2)=nit(2)+1
      mu = ddot(n1,x1,1,z1,1)+ddot(n2,x2,1,z2,1)+ddot(n1,s,1,w,1)
      g = mu + deltap*ddot(n1,dx1,1,z1,1)+ deltad*ddot(n1,dz1,1,x1,1)+ d
     *eltap*deltad*ddot(n1,dz1,1,dx1,1)+ deltap*ddot(n2,dx2,1,z2,1)+ del
     *tad*ddot(n2,dz2,1,x2,1)+ deltap*deltad*ddot(n2,dz2,1,dx2,1)+ delta
     *p*ddot(n1,ds,1,w,1)+ deltad*ddot(n1,dw,1,s,1) + deltap*deltad*ddot
     *(n1,ds,1,dw,1)
      mu = mu * ((g/mu)**3) /(dfloat(2*n1)+dfloat(n2))
      do23036 i=1,n1
      dsdw = ds(i)*dw(i)
      dr1(i)=d1(i)*(mu*(one/s(i)-one/x1(i))+ dx1(i)*dz1(i)/x1(i)-dsdw/s(
     *i))
23036 continue
23037 continue
      do23038 i=1,n2
      dr2(i)=d2(i)*(dx2(i)*dz2(i)/x2(i)-mu/x2(i))
23038 continue
23039 continue
      call dswap(p,rhs,1,dy,1)
      call dgemv('N',p,n1,one,a1,p,dr1,1,one,dy,1)
      call dgemv('N',p,n2,one,a2,p,dr2,1,one,dy,1)
      call dpotrs('U',p,1,ada,p,dy,p,info)
      call dgemv('T',p,n1,one,a1,p,dy,1,zero,u,1)
      deltap=big
      deltad=big
      do23040 i=1,n1
      dsdw = ds(i)*dw(i)
      dxdz1 = dx1(i)*dz1(i)
      dx1(i) = d1(i)*(u(i)-z1(i)+w(i))-dr1(i)
      ds(i) = -dx1(i)
      dz1(i) = -z1(i)+(mu - z1(i)*dx1(i) - dxdz1)/x1(i)
      dw(i) = -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
      if(dx1(i).lt.0)then
      deltap=dmin1(deltap,-x1(i)/dx1(i))
      endif
      if(ds(i).lt.0)then
      deltap=dmin1(deltap,-s(i)/ds(i))
      endif
      if(dz1(i).lt.0)then
      deltad=dmin1(deltad,-z1(i)/dz1(i))
      endif
      if(dw(i).lt.0)then
      deltad=dmin1(deltad,-w(i)/dw(i))
      endif
23040 continue
23041 continue
      call dgemv('T',p,n2,one,a2,p,dy,1,zero,u,1)
      do23050 i=1,n2
      dxdz2 = dx2(i)*dz2(i)
      dx2(i) = d2(i)*(u(i)-r2(i))-dr2(i)
      dz2(i) = -z2(i)+(mu - z2(i)*dx2(i) - dxdz2)/x2(i)
      if(dx2(i).lt.0)then
      deltap=dmin1(deltap,-x2(i)/dx2(i))
      endif
      if(dz2(i).lt.0)then
      deltad=dmin1(deltad,-z2(i)/dz2(i))
      endif
23050 continue
23051 continue
      deltap=dmin1(beta*deltap,one)
      deltad=dmin1(beta*deltad,one)
      endif
      call daxpy(n1,deltap,dx1,1,x1,1)
      call daxpy(n2,deltap,dx2,1,x2,1)
      call daxpy(n1,deltap,ds,1,s,1)
      call daxpy(p,deltad,dy,1,y,1)
      call daxpy(n1,deltad,dz1,1,z1,1)
      call daxpy(n2,deltad,dz2,1,z2,1)
      call daxpy(n1,deltad,dw,1,w,1)
      gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
      goto 23010
      endif
23011 continue
      call daxpy(n1,mone,w,1,z1,1)
      call dswap(n1,z1,1,x1,1)
      return
      end
      subroutine stepy2(n1,n2,p,a1,d1,a2,d2,b,ada,info)
      integer n1,n2,p,i,j,k,info
      double precision a1(p,n1),a2(p,n2),b(p),d1(n1),d2(n2),ada(p,p),zer
     *o
      parameter(zero = 0.0d0)
      do23056 j=1,p
      do23058 k=1,p
      ada(j,k)=zero
23058 continue
23059 continue
23056 continue
23057 continue
      do23060 i=1,n1
      call dsyr('U',p,d1(i),a1(1,i),1,ada,p)
23060 continue
23061 continue
      do23062 i=1,n2
      call dsyr('U',p,d2(i),a2(1,i),1,ada,p)
23062 continue
23063 continue
      call dposv('U',p,1,ada,p,b,p,info)
      return
      end
