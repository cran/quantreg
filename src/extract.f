c
c Extract: Subroutine to extract the non-diagonal structure and
c   entries from A stored in CSR format
c 
      subroutine extract(d,jd,id,dsub,jdsub,m,nnzd,nnzds,ierr)
      integer jd(nnzd),jdsub(nnzds),id(*),m,mp1,ierr,nnzd,nnzds
      double precision d(nnzd),dsub(nnzds)
c
c Call csrmsr in SPARSKIT2 to transform the storage format in d
c   from csr to msr
c
      call csrmsr(m,d,jd,id,dsub,jdsub,dsub,jdsub,nnzds,ierr)
      mp1 = m + 1
      do i=1,mp1
         jdsub(i) = jdsub(i) - mp1
      enddo
      return
      end
