c----------------------------------------------------------------------c
       subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A = row dimension of C
c ncol  = integer. The column dimension of B = column dimension of C
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c           
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in B.
c Note: 
c-------
c   The row dimension of B is not needed. However there is no checking 
c   on the condition that ncol(A) = nrow(B). 
c
c----------------------------------------------------------------------- 
      real*8 scal 
      logical values
      values = (job .ne. 0) 
      len = 0
      ic(1) = 1 
      ierr = 0
c     initialize array iw.
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c
      do 500 ii=1, nrow 
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
	    if (values) scal = a(ka)
	    jj   = ja(ka)
	    do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
 100	    continue
 200     continue
         do 201 k=ic(ii), len
	    iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
c-------------end-of-amub-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(*) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = A * Diag  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     scale each element 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c-----------------------------------------------------------------------
c-----------end-of-amudiag----------------------------------------------
      end 
c----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      real*8 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      do 500 ii=1, nrow
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
c     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
	    iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
c------------end of aplb ----------------------------------------------- 
c-----------------------------------------------------------------------
      end
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk,nnzao,ierr)
      real*8 a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1),nnzao,ierr
c----------------------------------------------------------------------- 
c Compressed Sparse Row   to      Modified - Sparse Row 
c                                 Sparse row with separate main diagonal
c----------------------------------------------------------------------- 
c converts a general sparse matrix a, ja, ia into 
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form 
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c----------------------------------------------------------------------- 
c
c on entry :
c---------
c a, ja, ia = matrix in csr format. note that the 
c	     algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c nnzao = the number of non-zero entries in ao, jao
c	 
c on return :
c-----------
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c	   +  ao(1:n) contains the diagonal of the matrix. 
c	   +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c	   +  jao(n+2:nnz) : their column indices
c	   +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+2:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i 
c	      in arrays ao, jao.
c	       here nnz = number of nonzero elements+1 
c ierr:
c    = -1 : length of ao, jao < itpr
c
c work arrays:
c------------
c wk	= real work array of length n
c iwk   = integer work array of length n+1
c
c notes: 
c------- 
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c--------
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c Modified by Pin Ng on June 11, 2002 to provide warning when 
c    iptr > nnzao+1
c-----------------------------------------------------------------------
      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1 
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
c     
c compute total length
c     
      iptr = n + ia(n+1) - icount
      if (iptr .gt. nnzao+1) then
         ierr = -1
         return
      endif
c     
c     copy backwards (to avoid collisions)
c     
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr-1
            endif
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i) 
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return	
c------------ end of subroutine csrmsr ---------------------------------
c----------------------------------------------------------------------- 
      end
