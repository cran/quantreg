c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c Subroutine to perform Cholesky factorization
c
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine chlfct(m,xlindx,lindx,invp,perm,iwork,nnzdsub,jdsub,
     &                  colcnt,nsuper,snode,xsuper,nnzlmax,nsubmax,
     &                  xlnz,lnz,id,jd,d,cachsz,tmpmax,level,
     &                  tmpvec,split,ierr,it,timewd)
      integer m,ierr,nsub,iwsiz,nnzdsub,nnzl,nsuper,nnzlmax,nsubmax,
     &        tmpsiz,cachsz,tmpmax,level,it,
     &        xlindx(*),lindx(*),invp(*),perm(*),iwork(*),jdsub(*),
     &        colcnt(*),snode(*),xsuper(*),xlnz(*),id(*),jd(*),split(*)
      double precision d(*),lnz(*),tmpvec(*),timewd(7)
      real gtimer,timbeg,timend
      external smxpy1,smxpy2,smxpy4,smxpy8
      external mmpy1,mmpy2,mmpy4,mmpy8
c
c Save the matrix structure from jdsub(m+2:nzzd+1),jdsub(1:m+1)
c   to lindx and xlindx because the matrix structure is destroyed by the 
c   minimum degree ordering routine
c
      do i = 1,m+1
         xlindx(i) = jdsub(i)
      enddo
      do i = 1,nnzdsub
         lindx(i) = jdsub(m+1+i)
      enddo
c
c Reorder the matrix using minimum degree ordering routine
c
      iwsiz = 4*m
      if(it .le. 1) then
         timbeg = gtimer()
         call ordmmd(m,xlindx,lindx,invp,perm,iwsiz,iwork,nsub,ierr)
         if (ierr .eq. -1) then
            ierr = 3
            go to 100
         endif
         timend = gtimer()
         timewd(1) = timewd(1) + timend - timbeg
c
c Call sfinit: Symbolic factorization initialization
c   to compute supernode partition and storage requirements
c   for symbolic factorization. New ordering is a postordering 
c   of the nodal elimination tree
c
         iwsiz = 7 * m + 3
         timbeg = gtimer()
         call sfinit(m,nnzdsub,jdsub(1),jdsub(m+2),perm,
     &            invp,colcnt,nnzl,nsub,nsuper,snode,xsuper,iwsiz,
     &            iwork,ierr)
         if (ierr .eq. -1) then
            ierr = 4
            go to 100
         endif
         if (nnzl .gt. nnzlmax) then
            ierr = 5
            go to 100
         endif
         if (nsub .gt. nsubmax) then
            ierr = 6
            go to 100
         endif
         timend = gtimer()
         timewd(2) = timewd(2) + timend - timbeg
      endif
c
c Call symfct: Perform supernodal symbolic factorization
c
c      iwsiz = nsuper + 2 * m + 1
      timbeg = gtimer()
      call symfct(m,nnzdsub,jdsub(1),jdsub(m+2),perm,invp,
     &            colcnt,nsuper,xsuper,snode,nsub,xlindx,lindx,
     &            xlnz,iwsiz,iwork,ierr)
      if (ierr .eq. -1) then
         ierr = 7
         go to 100
      endif
      if (ierr .eq. -2) then
         ierr = 8
         go to 100
      endif
      timend = gtimer()
      timewd(3) = timewd(3) + timend - timbeg
c
c Call inpnv: Input numerical values into data structures of L
c
      iwsiz = m
      timbeg = gtimer()
      call inpnv(m,id,jd,d,perm,invp,nsuper,xsuper,xlindx,lindx,
     &           xlnz,lnz,iwork)
      timend = gtimer()
      timewd(4) = timewd(4) + timend - timbeg
c
c Call bfinit: Initialization for block factorization
c
      timbeg = gtimer()
      call bfinit(m,nsuper,xsuper,snode,xlindx,lindx,cachsz,tmpsiz,
     &            split)
      if (tmpsiz .gt. tmpmax) then
         ierr = 9
         go to 100
      endif
      timend = gtimer()
      timewd(5) = timewd(5) + timend - timbeg
c
c Call blkfct: Numerical factorization
c
      iwsiz = 2 * m + 2 * nsuper
      timbeg = gtimer()
      if (level .eq. 1) then
         call blkfct(m,nsuper,xsuper,snode,split,xlindx,lindx,xlnz,
     &               lnz,iwsiz,iwork,tmpsiz,tmpvec,ierr,mmpy1,smxpy1)
      elseif (level .eq. 2) then
         call blkfct(m,nsuper,xsuper,snode,split,xlindx,lindx,xlnz,
     &               lnz,iwsiz,iwork,tmpsiz,tmpvec,ierr,mmpy2,smxpy2)
      elseif (level .eq. 4) then
         call blkfct(m,nsuper,xsuper,snode,split,xlindx,lindx,xlnz,
     &               lnz,iwsiz,iwork,tmpsiz,tmpvec,ierr,mmpy4,smxpy4)
      elseif (level .eq. 8) then
         call blkfct(m,nsuper,xsuper,snode,split,xlindx,lindx,xlnz,
     &               lnz,iwsiz,iwork,tmpsiz,tmpvec,ierr,mmpy8,smxpy8)
      endif
      if (ierr .eq. -1) then
         ierr = 10
         go to 100
      elseif (ierr .eq. -2) then
         ierr = 11
         go to 100
      elseif (ierr .eq. -3) then
         ierr = 12
         go to 100
      endif
  100 continue
      timend = gtimer()
      timewd(6) = timewd(6) + timend - timbeg
      return
      end
c
