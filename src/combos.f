C Output from Public domain Ratfor, version 1.0
      subroutine combin(r,n,m,a,c,e,last)
      integer r,n,m,t,k,j,m0,mj
      integer a(n,m),c(r),e(r),last(r)
      logical odd
      m0 = r-n
      t = n+1
      k = 1
      j = 0
      c(1) = 0
23000 continue
      j = j + 1
      c(j) = j
      e(j) = j - 1
      if(odd(j))then
      last(j) = m0 + j
      else
      last(j) = j + 1
      endif
      if(j .eq. n)then
      goto 23002
      endif
23001 goto 23000
23002 continue
      do23007 i = 1,n
      a(i,1) = c(i)
23007 continue
23008 continue
      if(n .lt. r)then
23011 continue
      k = k + 1
      s = c(j)
      mj = m0 + j
      e(n+1) = n
      if(odd(j))then
      if(c(j) .eq. mj)then
      c(j) = c(j-1) + 1
      last(j+1) = c(j) + 1
      else
      c(j) = c(j) + 1
      endif
      else
      if(c(j) .eq. c(j-1) + 1)then
      c(j) = mj
      else
      last(j+1) = c(j)
      c(j) = c(j) - 1
      endif
      endif
      if(c(j) .eq. last(j))then
      last(j) = s
      e(j+1) = e(j)
      e(j) = j-1
      endif
      if( (j .lt. n) .and. (c(j) .eq. mj))then
      t = j
      j = e(t+1)
      e(t+1) = t
      else
      if(t .eq. j)then
      t = t + 1
      endif
      if(t .lt. e(n+1))then
      j = t
      else
      j = e(n+1)
      endif
      endif
      do23028 i = 1,n
      a(i,k) = c(i)
23028 continue
23029 continue
      if(j .eq. 0)then
      goto 23013
      endif
23012 goto 23011
23013 continue
      endif
      return
      end
      logical function odd(j)
      integer j
      odd = (mod(j,2) .eq. 1)
      return
      end
