c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine srqfnc(n1,m,nnza1, a1,ja1,ia1, ao1,jao1,iao1, n2,nnza2, ! 11
     &                  a2,ja2,ia2, ao2,jao2,iao2, nnzdmax, d,jd,id,     ! 21
     &                  dsub,jdsub, nnzemax, e,je,ie, nnzgmax, g,jg,ig,  ! 31
     &                  nnzhmax, h,jh,ih, nsubmax, lindx,xlindx,nnzlmax, ! 39
     &                  lnz,xlnz, iw,iwmax,iwork, xsuper, tmpmax,tmpvec, ! 47
     &                  maxn1n2, ww1,wwm,wwn1,wwn2, cachsz, level,x1,x2, ! 56
     &                  s,u,c1,c2, sm_tn_Lrg,   y, ierr,maxit, timewd)   ! 65 ( = MAX_ARGS !)
      integer nnza1,nnza2,m,n1,n2,nnzdmax,nnzemax,nnzgmax,nnzhmax,iwmax,
     &        nnzlmax,nsubmax,cachsz,level,tmpmax,ierr,maxit,maxn1n2,
     &        ja1(nnza1),jao1(nnza1),ja2(nnza2),jao2(nnza2),
     &        jdsub(nnzhmax+1),jd(nnzdmax),ia1(n1+1),iao1(m+1),
     &        ia2(n2+1),iao2(m+1),id(m+1),lindx(nsubmax),xlindx(m+1),
     &        iw(m,5),xlnz(m+1),iwork(iwmax),xsuper(m+1),je(nnzemax),
     &        ie(m+1),jg(nnzgmax),ig(m+1),jh(nnzhmax),ih(m+1)
      double precision sm_tn_Lrg(3), !  := c(small, tiny, Large)
     &                 a1(nnza1),ao1(nnza1),a2(nnza2),ao2(nnza2),
     &                 dsub(nnzhmax+1),d(nnzdmax),g(nnzgmax),
     &                 h(nnzhmax),lnz(nnzlmax),c1(n1),c2(n2),y(m),
     &                 ww1(maxn1n2),wwm(m,6),tmpvec(tmpmax),
     &                 wwn1(n1,10),wwn2(n2,7),x1(n1),x2(n2),s(n1),
     &                 u(n1),e(nnzemax)
      double precision timewd(7)
      parameter (beta=9.995d-1, one=1.0d0, zero=0.0d0)
      call slpfnc(n1,m,nnza1,a1,ja1,ia1,ao1,jao1,iao1,n2,nnza2,
     &          a2,ja2,ia2,ao2,jao2,iao2,nnzdmax,d,jd,id,dsub,
     &          jdsub,nsubmax,lindx,xlindx,nnzlmax,lnz,xlnz,iw(1,1),
     &          iw(1,2),iwmax,iwork,iw(1,3),iw(1,4),xsuper,iw(1,5),
     &          tmpmax,tmpvec,wwm(1,2),wwm(1,3),cachsz,level,x1,x2,s,u,
     &          c1,c2,y,wwm(1,1),wwn2(1,1),wwn1(1,1),
     &          wwn2(1,2),wwn1(1,2),wwn1(1,3),wwn2(1,3),nnzemax,e,je,
     &          ie,nnzgmax,g,jg,ig,nnzhmax,h,jh,ih,wwm(1,4),wwn1(1,4),
     &          wwn2(1,4),wwn1(1,5),wwn1(1,6),wwn2(1,5),wwn1(1,7),
     &          wwn1(1,8),wwn2(1,6),wwn1(1,9),wwn1(1,10),wwn2(1,7),
     &     maxn1n2,ww1,wwm(1,5),wwm(1,6), sm_tn_Lrg(1), ierr,maxit,
     &     timewd, sm_tn_Lrg(2), sm_tn_Lrg(3))
      return
      end

c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine slpfnc(n1,m,nnza1,a1,ja1,ia1,ao1,jao1,iao1,n2,nnza2,
     &                a2,ja2,ia2,ao2,jao2,iao2,nnzdmax,d,jd,id,dsub,
     &                jdsub,nsubmax,lindx,xlindx,nnzlmax,lnz,xlnz,invp,
     &                perm,iwmax,iwork,colcnt,snode,xsuper,split,
     &                tmpmax,tmpvec,rhs,newrhs,cachsz,level,x1,x2,s,u,
     &                c1,c2,y,b,r2,z1,
     &                z2,w,q1,q2,nnzemax,e,je,
     &                ie,nnzgmax,g,jg,ig,nnzhmax,h,jh,ih,dy,dx1,
     &                dx2,ds,dz1,dz2,dw,
     &                dxdz1,dxdz2,dsdw,xi1,xi2,
     &                maxn1n2, ww1,ww2,ww3, small, ierr,maxit, timewd,
     &                tiny, Large)
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c Sparse implentation of LMS's interior point method via
c    Ng-Peyton's sparse Cholesky factorization for sparse
c    symmetric positive definite
c INPUT:
c     n1 -- the number of row in the coefficient matrix A1'
c     m -- the number of column in the coefficient matrix A1'
c     nnza1 -- the number of non-zero elements in A'
c     a1 -- an nnza1-vector of non-zero values of the design
c          matrix (A1') stored in csr format
c     ja1 -- an nnza1-vector of indices of the non-zero elements of
c           the coefficient matrix
c     ia1 -- an (n1+1)-vector of pointers to the begining of each
c           row in a1 and ja1
c     ao1 -- an nnza1-vector of work space for the transpose of
c           the design matrix stored in csr format or the
c           design matrix stored in csc format
c     jao1 -- an nnza1-vector of work space for the indices of the
c            transpose of the design matrix
c     iao1 -- an (n1+1)-vector of pointers to the begining of each
c            column in ao1 and jao1
c     n2 -- the number of row in the constraint matrix A2'
c     nnza2 -- the number of non-zero elements in A2'
c     a2 -- an nnza2-vector of non-zero values of the contraint
c          matrix (A2') stored in csr format
c     ja2 -- an nnza2-vector of indices of the non-zero elements of
c           the constraint matrix
c     ia2 -- an (n2+1)-vector of pointers to the begining of each
c           row in a2 and ja2
c     ao2 -- an nnza2-vector of work space for the transpose of
c           the constraint matrix stored in csr format or the
c           constraint matrix stored in csc format
c     jao2 -- an nnza2-vector of work space for the indices of the
c            transpose of the constraint matrix
c     iao2 -- an (n2+1)-vector of pointers to the begining of each
c            column in ao2 and jao2
c     nnzdmax -- upper bound of the non-zero elements in A1A1'
c     d -- an nnzdmax-vector of non-zero values used to store
c          the transpose of the design matrix multiplied by the design
c          matrix (A1A1') stored in csr format;
c          also used to store A1Q1^(-1) and A2Q2^(-1) later
c     jd -- an nnzdmax-vector of indices in d
c     id -- an (m+1)-vector of pointers to the begining of each
c           row in d and jd
c     dsub -- the values of d excluding the diagonal elements
c     jdsub -- the indices to dsub
c     nsubmax -- upper bound of the dimension of lindx
c     lindx -- an nsub-vector of interger which contains, in
c           column major order, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an (m+1)-vector of integer of pointers for lindx
c     nnzlmax -- the upper bound of the non-zero entries in
c                L stored in lnz, including the diagonal entries
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     invp -- an n1-vector of integer of inverse permutation
c             vector
c     perm -- an n1-vector of integer of permutation vector
c     iw -- integer work array of length m
c     iwmax -- upper bound of the general purpose integer
c              working storage iwork; set at 7*m+3
c     iwork -- an iwsiz-vector of integer as work space
c     colcnt -- array of length m, containing the number of
c               non-zeros in each column of the factor, including
c               the diagonal entries
c     snode -- array of length m for recording supernode
c              membership
c     xsuper -- array of length m+1 containing the supernode
c               partitioning
c     split -- an m-vector with splitting of supernodes so that
c              they fit into cache
c     tmpmax -- upper bound of the dimension of tmpvec
c     tmpvec -- a tmpmax-vector of temporary vector
c     rhs -- m-vector to store the rhs
c     newrhs -- extra work vector for right-hand side and
c               solution
c     cachsz -- size of the cache (in kilobytes) on the target
c               machine
c     level -- level of loop unrolling while performing numerical
c              factorization
c     x1 -- an n1-vector, the initial feasible solution for the primal
c           solution that corresponds to the design matrix A1'
c     x2 -- an n2-vector, the initial feasible solution for the primal
c           solution that corresponds to the constraint matrix A2'
c     s -- an n1-vector
c     u -- an n1-vector of the upper bound for x1
c     c1 -- an n1-vector in the primal; negative response in the
c           regression quantile setting
c     c2 -- an n2-vector, the negative rhs of the inequality constraint
c     y -- an m-vector, the initial dual solution
c     b -- an n1-vector, usualy the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c     r2 -- an n2-vector of residuals
c     z1 -- an n1-vector of the dual slack variable
c     z2 -- an n2-vector
c     w -- an n-vector
c     q1 -- an n1-vector of work array containing the diagonal
c          elements of the Q1^(-1) matrix
c     q2 -- an n2-vector of work array containing the diagonal
c          elements of the Q2^(-1) matrix
c     e -- an nnzdmax-vector containing the non-zero entries of
c          A1Q1^(-1)A1' stored in csr format
c     je -- an nnzdmax-vector of indices for e
c     ie -- an (m+1)-vector of pointers to the begining of each
c           row in e and je
c     nnzgmax -- upper bound of the non-zero elements in g,jg
c     g -- an nnzgmax-vector containing the non-zero entries of
c          A2Q2^(-1)A2' stored in csr format
c     jg -- an nnzgmax-vector of indices for g
c     ig -- an (m+1)-vector of pointers to the begining of each
c           row in g and jg
c     nnzhmax -- upper bound of the non-zero elements in h,jh
c     h -- an nnzhmax-vector containing the non-zero entries of
c          AQ^(-1)A' stored in csr format
c     jh -- an nnzhmax-vector of indices for h
c     ih -- an (m+1)-vector of pointers to the begining of each
c           row in h and jh
c     dy -- an m-vector of work array
c     dx1 -- an n1-vector of work array
c     dx2 -- an n2-vector of work array
c     ds -- an n1-vector of work array
c     dz1 -- an n1-vector of work array
c     dz2 -- an n2-vector of work array
c     dw -- an n1-vector of work array
c     dxdz1 -- an n1-vector of work array
c     dxdz2 -- an n2-vector of work array
c     dsdw -- an n1-vector of work arry
c     xi1 -- an n1-vector of work array
c     xi2 -- an n2-vector of work array
c     xinv1 -- an n1-vector of work array
c     xinv2 -- an n2-vector of work array
c     sinv -- work array
c     maxn1n2 -- max(n1,n2)
c     ww1 -- an maxn1n2-vector of work array
c     ww2 -- an m-vector of work array
c     ww3 -- an m-vector of work array
c     small -- convergence tolerance for inetrior algorithm
c     ierr -- error flag
c       1 -- insufficient storage when calling extract;
c       3 -- insufficient storage in iwork when calling ordmmd;
c       4 -- insufficient storage in iwork when calling sfinit;
c       5 -- nnzl > nnzlmax when calling sfinit
c       6 -- nsub > nsubmax when calling sfinit
c       7 -- insufficient work space in iwork when calling symfct
c       8 -- inconsistancy in input when calling symfct
c       9 -- tmpsiz > tmpmax when calling symfct; increase tmpmax
c       10 -- nonpositive diagonal encountered when calling
c            blkfct
c       11 -- insufficient work storage in tmpvec when calling
c            blkfct
c       12 -- insufficient work storage in iwork when calling
c            blkfct
c       13 -- nnzd > nnzdmax in e,je when calling amub
c       14 -- nnzd > nnzdmax in g,jg when calling amub
c       15 -- nnzd > nnzdmax in h,jh when calling aplb
c     maxit -- upper limit of the iteration; on return holds the
c              number of iterations
c     timewd -- amount of time to execute this subroutine
c     tiny  -- tiny number; values below tiny * max(diag) are replaced by 'Large';
c               was 10^{-30} hardcoded
c     Large -- Large number ("Infinite") to replace tiny diagonal entries in Cholesky;
c               was 10^{128}
c OUTPUT:
c     y -- an m-vector of primal solution
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer nnza1,nnza2,m,n1,n2,maxn1n2,nsuper,nnzdmax,nnzemax,
     &        nnzgmax,
     &        nnzhmax,iwmax,nnzlmax,nsubmax,cachsz,level,tmpmax,ierr,
     &        maxit,it,nnzdsub,nnzd,
     &        ja1(nnza1),jao1(nnza1),ja2(nnza2),jao2(nnza2),
     &        jdsub(nnzhmax+1),jd(nnzdmax),ia1(n1+1),iao1(m+1),
     &        ia2(n2+1),iao2(m+1),id(m+1),lindx(nsubmax),xlindx(m+1),
     &        invp(m),perm(m),xlnz(m+1),iwork(iwmax),colcnt(m),
     &        snode(m),xsuper(m+1),split(m),je(nnzemax),ie(m+1),
     &        jg(nnzgmax),ig(m+1),jh(nnzhmax),ih(m+1)
      double precision ddot,gap,zero,one,beta,small,deltap,deltad,mu,g1,
     &                 a1(nnza1),ao1(nnza1),a2(nnza2),ao2(nnza2),
     &                 dsub(nnzhmax+1),d(nnzdmax),lnz(nnzlmax),c1(n1),
     &                 c2(n2),b(m),rhs(m),newrhs(m),y(m),tmpvec(tmpmax),
     &                 r2(n2),z1(n1),z2(n2),w(n1),x1(n1),
     &                 x2(n2),s(n1),u(n1),q1(n1),q2(n2),e(nnzemax),
     &                 g(nnzgmax),h(nnzhmax),dy(m),dx1(n1),dx2(n2),
     &                 ds(n1),dz1(n1),dz2(n2),dw(n1),dxdz1(n1),
     &                 dxdz2(n2),dsdw(n1),
     &                 xi1(n1),xi2(n2),ww1(maxn1n2),ww2(m),ww3(m)
      double precision timewd(7), tiny,Large
      real gtimer,timbeg,timend
      external smxpy1,smxpy2,smxpy4,smxpy8
      external mmpy1,mmpy2,mmpy4,mmpy8
      parameter (beta=9.995d-1, one=1.0d0, zero=0.0d0)
      do i = 1,7
         timewd(i) = 0.0
      enddo
      it = 0
c
c  Compute the initial gap
c
       gap = ddot(n1,z1,1,x1,1) + ddot(n2,z2,1,x2,1) + ddot(n1,w,1,s,1)
c
c  Start iteration
c
   20 continue
      if(gap .lt. small .or. it .gt. maxit) goto 30
         it = it + 1
c
c  Create the diagonal matrix Q1^(-1) stored in q1 as an n1-vector,
c  the diagonal matrix Q2^(-1) stored in q2 as an n2-vector,
c  and store the residuals in r1 in ds, and r3 in dy temporarily,
c  and r2 in r2 permanently
c
c    Call amux to obtain A1x1 and store the value in ww2
c
         call amux(m,x1,ww2,ao1,jao1,iao1)
c
c    Call amux to obtain A2x2 and store the value in ww3
c
         call amux(m,x2,ww3,ao2,jao2,iao2)
c
c    Store A2'y temporarily in r2
c
         call amux(n2,y,r2,a2,ja2,ia2)
         do i=1,n1
            q1(i) = one/(z1(i)/x1(i)+w(i)/s(i))
            ds(i) = z1(i) - w(i)
         enddo
         do i=1,n2
            q2(i) = x2(i)/z2(i)
            r2(i) = c2(i)-r2(i)
         enddo
         do i=1,m
            dy(i) = b(i) - ww2(i) - ww3(i)
         enddo
c
c  Obtain AQA = A1Q1^(-1)A1' + A2Q2^(-1)A2' in 5 steps
c
c  Step1: Obtain A1Q1^(-1) and store the values in d,jd,id in csr format
c         Also compute A1Q1^(-1)r1 and store the values in ww2 to be used
c         to generate r3;
c  Step2: Compute A1Q1^(-1)A1' and store the values in e,je,ie
c  Step3: Obtain A2Q2^(-1) and store the values in d,jd,id in csr format
c         Also compute A2Q2^(-1)r2 and store the values in in ww3 to
c         be used to generate r3;
c  Step4: Compute A2Q2^(-1)A2' and store the value in g,jg,ig
c  Step5: Compute AQA and store the values in h,jh,ih
c
c    Step 1
c
         call amudia(m,1,ao1,jao1,iao1,q1,d,jd,id)
         call amux(m,ds,ww2,d,jd,id)
c
c   Step 2
c
         call amub(m,m,1,d,jd,id,a1,ja1,ia1,e,je,ie,nnzemax,iwork,ierr)
         if (ierr .ne. 0) then
            ierr = 13
            go to 100
         endif
c
c   Step 3
c
         call amudia(m,1,ao2,jao2,iao2,q2,d,jd,id)
         call amux(m,r2,ww3,d,jd,id)
c
c   Step 4
c
         call amub(m,m,1,d,jd,id,a2,ja2,ia2,g,jg,ig,nnzgmax,iwork,ierr)
         if (ierr .ne. 0) then
            ierr = 14
            go to 100
         endif
c
c   Step 5
c
         call aplb(m,m,1,e,je,ie,g,jg,ig,h,jh,ih,nnzhmax,iwork,ierr)
         if (ierr .ne. 0) then
            ierr = 15
            go to 100
         endif
c
c  Generate rhs = r3 + A1Q1^(-1)r1 + A2Q2^(-1)r2 and store in rhs
c
         do i = 1,m
            rhs(i) = dy(i) + ww2(i) + ww3(i)
         enddo
c
c  Extract the non-diagonal structure of h,jh,ih and store in dsub,jdsub
c
         nnzd = ih(m+1) - 1
         nnzdsub = nnzd - m
         call extract(h,jh,ih,dsub,jdsub,m,nnzhmax,nnzhmax+1,ierr)
         if (ierr .ne. 0) then
            ierr = 1
            go to 100
         endif
c
c  Compute dy = (AQ^(-1)A')^(-1)rhs; result returned via dy
c
c Call chlfct to perform Cholesky's decomposition of h,jh,ih
c
         call chlfct(m,xlindx,lindx,invp,perm,iwork,nnzdsub,jdsub,
     &            colcnt,nsuper,snode,xsuper,nnzlmax,nsubmax,xlnz,lnz,
     &            ih,jh,h,cachsz,tmpmax,level,tmpvec,split,ierr,it,
     &            timewd, tiny,Large)

         if (ierr .ne. 0) go to 100
c
c Call blkslv: Numerical solution for the new rhs stored in rhs
c
         do i = 1,m
            newrhs(i) = rhs(perm(i))
         enddo
         timbeg = gtimer()
         call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
         timend = gtimer()
         timewd(7) = timewd(7) + timend - timbeg
         do i = 1,m
            dy(i) = newrhs(invp(i))
         enddo
c
c  Compute dx1 = Q1^(-1)(A1'dy - r1), ds = -dx1, dz1, dz2  and dw
c
         call amux(n1,dy,dx1,a1,ja1,ia1)
         call amux(n2,dy,dx2,a2,ja2,ia2)
         do i = 1,n1
            dx1(i) = q1(i) * (dx1(i) - ds(i))
            ds(i) = -dx1(i)
            dz1(i) = -z1(i) * (one + dx1(i) / x1(i))
            dw(i) = -w(i) * (one + ds(i) / s(i))
         enddo
         do i = 1,n2
            dx2(i) = q2(i) * (dx2(i) - r2(i))
            dz2(i) = -z2(i) * (one + dx2(i) / x2(i))
         enddo
c
c Compute the maximum allowable step lengths
c
         call boundc(x1,dx1,x2,dx2,s,ds,z1,dz1,z2,dz2,w,dw,n1,n2,
     &               beta,deltap,deltad)
         if (deltap * deltad .lt. one) then
c
c Update mu
c
            mu = ddot(n1,z1,1,x1,1) + ddot(n2,z2,1,x2,1)
     &           + ddot(n1,w,1,s,1)
            g1 = mu + deltap*ddot(n1,z1,1,dx1,1)
     &          + deltad*ddot(n1,dz1,1,x1,1)
     &          + deltad*deltap*ddot(n1,dz1,1,dx1,1)
     &          + deltap*ddot(n2,z2,1,dx2,1)
     &          + deltad*ddot(n2,dz2,1,x2,1)
     &          + deltad*deltap*ddot(n2,dz2,1,dx2,1)
     &          + deltap*ddot(n1,w,1,ds,1)
     &          + deltad*ddot(n1,dw,1,s,1)
     &          + deltad*deltap*ddot(n1,dw,1,ds,1)
            mu = mu*((g1/mu)**3)/(2.d0*dble(n1)+dble(n2))
c
c Compute dx1dz1, dx2dz2 and dsdw
c
            do i = 1,n1
                dxdz1(i) = dx1(i)*dz1(i)
                dsdw(i) = ds(i)*dw(i)
                xi1(i) = dxdz1(i)/x1(i) - dsdw(i)/s(i)
     &                  - mu * (one/x1(i) - one/s(i))

                ww1(i) = q1(i) * xi1(i)
            enddo
c
c Compute A1Q1^(-1)(X1^(-1)*dx1dz1 - S^(-1)*dsdw - mu(X1^(-1) - S^(-1))) and
c store it in ww2 temporarily
c
            call amux(m,ww1,ww2,ao1,jao1,iao1)
            do i = 1,n2
                dxdz2(i) = dx2(i)*dz2(i)
                xi2(i) = (dxdz2(i) - mu)/x2(i)
                ww1(i) = q2(i) * xi2(i)
            enddo
c
c Compute A2Q2^(-1)(X2^(-1)*dx2dz2 - mu X2^(-1)) and store it in ww3
c temporarily
c
            call amux(m,ww1,ww3,ao2,jao2,iao2)
            do i = 1,m
               rhs(i) = rhs(i) + ww2(i) + ww3(i)
            enddo
c
c
c Compute (AQ^(-1)A')^(-1)rhs and return the result in dy
c
c Call blkslv: Numerical solution for the new rhs stored in rhs
c
            do i = 1,m
               newrhs(i) = rhs(perm(i))
            enddo
            timbeg = gtimer()
            call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
            timend = gtimer()
            timewd(7) = timewd(7) + timend - timbeg
            do i = 1,m
               dy(i) = newrhs(invp(i))
            enddo
c
c  Compute dx1=Q1^(-1)(A1'dy-X1^(-1)*dx1dz1-S^(-1)*dsdw
c  -mu*(X1^(-1)-S^(-1))-r1), ds = -dx1, dz1, dz2  and dw
c
            call amux(n1,dy,dx1,a1,ja1,ia1)
            call amux(n2,dy,dx2,a2,ja2,ia2)
            do i = 1,n1
               dx1(i) = q1(i) * (dx1(i) - xi1(i) - z1(i) + w(i))
               ds(i) = -dx1(i)
               dz1(i) = -z1(i) + (mu - z1(i)*dx1(i)
     &                  - dxdz1(i))/x1(i)
               dw(i) = -w(i) + (mu - w(i)*ds(i) - dsdw(i))/s(i)
            enddo
            do i = 1,n2
               dx2(i) = q2(i) * (dx2(i) - xi2(i) - r2(i))
               dz2(i) = -z2(i) + (mu - z2(i)*dx2(i) -
     &                  dxdz2(i))/x2(i)
            enddo
c
c Compute the maximum allowable step lengths
c
            call boundc(x1,dx1,x2,dx2,s,ds,z1,dz1,z2,dz2,w,dw,n1,n2,
     &                  beta,deltap,deltad)
         endif
c
c Take the step
c
         call daxpy(n1,deltap,dx1,1,x1,1)
         call daxpy(n2,deltap,dx2,1,x2,1)
         call daxpy(n1,deltap,ds,1,s,1)
         call daxpy(n1,deltad,dw,1,w,1)
         call daxpy(n1,deltad,dz1,1,z1,1)
         call daxpy(n2,deltad,dz2,1,z2,1)
         call daxpy(m,deltad,dy,1,y,1)
         gap = ddot(n1,z1,1,x1,1) + ddot(n2,z2,1,x2,1) +
     &         ddot(n1,w,1,s,1)
      goto 20
   30 continue
  100 continue
      maxit = it
      return
      end
