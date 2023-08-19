C Output from Public domain Ratfor, version 1.05
      subroutine qselect(n,x,q)
      integer n,k,l,r
      double precision x(n),q
      k=nint(q*n)
      l=1
      r=n
      call select(n,x,l,r,k)
      q=x(k)
      return
      end
      recursive subroutine select(n,x,l,r,k)
      integer n,m,l,r,k,ll,rr,i,j,mmax
      double precision x(n),z,s,d,t,fm,cs,cd
      parameter(cs = 0.5d0)
      parameter(cd = 0.5d0)
      parameter(mmax = 600)
23000 if(r.gt.l)then
      if(r-l.gt.mmax)then
      m=r-l+1
      i=k-l+1
      fm = dble(m)
      z=log(fm)
      s=cs*exp(2*z/3)
      d=cd*sqrt(z*s*(m-s)/fm)*sign(1,i-m/2)
      ll=max(l,nint(k-i*s/fm + d))
      rr=min(r,nint(k+(m-i)*s/fm + d))
      call select(n,x,ll,rr,k)
      endif
      t=x(k)
      i=l
      j=r
      call dswap(1,x(l),1,x(k),1)
      if(x(r).gt.t)then
      call dswap(1,x(r),1,x(l),1)
      endif
23006 if(i.lt.j)then
      call dswap(1,x(i),1,x(j),1)
      i=i+1
      j=j-1
23008 if(x(i).lt.t)then
      i=i+1
      goto 23008
      endif
23009 continue
23010 if(x(j).gt.t)then
      j=j-1
      goto 23010
      endif
23011 continue
      goto 23006
      endif
23007 continue
      if(x(l).eq.t)then
      call dswap(1,x(l),1,x(j),1)
      else
      j=j+1
      call dswap(1,x(j),1,x(r),1)
      endif
      if(j.le.k)then
      l=j+1
      endif
      if(k.le.j)then
      r=j-1
      endif
      goto 23000
      endif
23001 continue
      return
      end
