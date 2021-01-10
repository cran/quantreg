C Output from Public domain Ratfor, version 1.05
      subroutine grexp(n, x, a)
      integer i,n
      double precision x(n),a
      call fseedi()
      do23000 i = 1,n
      call frexp(x(i), a)
23000 continue
23001 continue
      call fseedo()
      return
      end
