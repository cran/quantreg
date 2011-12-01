C Output from Public domain Ratfor, version 1.0
      subroutine kuantile(k,m,n,x)
      integer i,j,k(m),m,n
      double precision x(n)
      j = 0
      do23000 i = 1,m
      call dsel05(k(i)-j,n-j,x(j+1))
      j = k(i)
23000 continue
23001 continue
      return
      end
