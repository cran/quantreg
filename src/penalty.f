C Output from Public domain Ratfor, version 1.0
      subroutine penalty(n,m,q,x,y,bnd,tlist,tlptr,tlend,rax,jax,ned,eps
     *,ierr)
      integer n,m,q,lp,lpl,ned,ierr
      integer bnd(n),tlist(q),tlptr(q),tlend(n),n4(4),p4(4),jax(m)
      double precision x(n),y(n),rax(m),eps
      double precision x4(4),y4(4),g4(4)
      logical orient
      ned = 0
      do23000 i=1,n
      lpl = tlend(i)
      lp = lpl
23002 continue
      lp = tlptr(lp)
      j = iabs(tlist(lp))
      if(j .gt. i)then
      n4(1) = i
      n4(2) = j
      call fadjs(n4,n,q,tlist,tlptr,tlend)
      if(bnd(i)*bnd(j) .eq. 0)then
      ned = ned + 1
      do23009 k = 1,4
      x4(k) = x(n4(k))
      y4(k) = y(n4(k))
23009 continue
23010 continue
      if(orient(x4,y4))then
      call iswap(1,n4(3),1,n4(4),1)
      call dswap(1,x4(3),1,x4(4),1)
      call dswap(1,y4(3),1,y4(4),1)
      endif
      call ggap(x4,y4,g4,eps,ierr)
      if(ierr .eq. 1)then
      return
      endif
      call srtpai(n4,1,p4,1,4)
      do23015 k = 1,4
      rax((ned - 1)*4 + k) = g4(p4(k))
      jax((ned - 1)*4 + k) = n4(p4(k))
23015 continue
23016 continue
      if(ned*4 .gt. m)then
      return
      endif
      endif
      endif
      if(lp .eq. lpl)then
      goto 23004
      endif
23003 goto 23002
23004 continue
23000 continue
23001 continue
      return
      end
      logical function orient(x,y)
      double precision x(4), y(4)
      orient = (y(2) -y(1))*(x(3)-x(4))+(x(1)-x(2))*(y(3)-y(4)) .gt. 0
      return
      end
      subroutine fadjs(n4,n,q,tlist,tlptr,tlend)
      integer n,q,vp,vpl,v,v0,match
      integer n4(4),tlist(q),tlptr(q),tlend(n)
      match = 0
      vpl = tlend(n4(1))
      vp = vpl
      k = 0
23021 continue
      k = k+1
      vp = tlptr(vp)
      v = tlist(vp)
      if(k.gt.1 .and. iabs(v) .eq. n4(2))then
      n4(3) = iabs(v0)
      match = 1
      goto 23022
      endif
      if(match .gt. 0)then
      n4(4) = iabs(v)
      goto 23023
      endif
      v0 = v
23022 goto 23021
23023 continue
      return
      end
      subroutine ggap(x,y,g,eps,ierr)
      double precision x(4),y(4),g(4),w(2,4),h(2),d1,d2,eps
      d1 = -x(2) * y(1) + x(3) * y(1) + x(1) * y(2) - x(3) * y(2) - x(1)
     * * y(3) + x(2) * y(3)
      d2 = -x(2) * y(1) + x(4) * y(1) + x(1) * y(2) - x(4) * y(2) - x(1)
     * * y(4) + x(2) * y(4)
      if(dabs(d1) .lt. eps .or. dabs(d2) .lt. eps)then
      ierr = 1
      return
      endif
      h(1) = -(y(1) - y(2))
      h(2) = (x(1) - x(2))
      w(1, 1) = (y(2) - y(3))/d1 - (y(2) - y(4))/d2
      w(2, 1) = (x(3) - x(2))/d1 - (x(4) - x(2))/d2
      w(1, 2) = (y(3) - y(1))/d1 - (y(4) - y(1))/d2
      w(2, 2) = (x(1) - x(3))/d1 - (x(1) - x(4))/d2
      w(1, 3) = (y(1) - y(2))/d1
      w(2, 3) = (x(2) - x(1))/d1
      w(1, 4) = (y(2) - y(1))/d2
      w(2, 4) = (x(1) - x(2))/d2
      do23030 i = 1,4
      g(i) = h(1)*w(1,i)+h(2)*w(2,i)
23030 continue
23031 continue
      ierr = 0
      return
      end
