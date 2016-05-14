c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine srqfn(n,m,nnza,a,ja,ia,ao,jao,iao,nnzdmax,d,jd,id,
     &                  dsub,jdsub,nnzemax,e,je,ie,nsubmax,lindx,xlindx,
     &                  nnzlmax,lnz,xlnz,iw,iwmax,iwork,xsuper,tmpmax,
     &                  tmpvec,wwm,wwn,cachsz,level,x,s,u,c,y,b,small,
     &                  ierr,maxit,timewd)
      integer nnza,m,n,nnzdmax,nnzemax,iwmax,
     &        nnzlmax,nsubmax,cachsz,level,tmpmax,ierr,maxit,
     &        ja(nnza),jao(nnza),jdsub(nnzemax+1),jd(nnzdmax),
     &        ia(n+1),iao(m+1),id(m+1),lindx(nsubmax),xlindx(m+1),
     &        iw(m,5),xlnz(m+1),iwork(iwmax),xsuper(m+1),je(nnzemax),
     &        ie(m+1)
      double precision small,
     &                 a(nnza),ao(nnza),dsub(nnzemax+1),d(nnzdmax),
     &                 lnz(nnzlmax),c(n),y(m),wwm(m,3),tmpvec(tmpmax),
     &                 wwn(n,14),x(n),s(n),u(n),e(nnzemax),b(m)
      double precision timewd(7)
      call slpfn(n,m,nnza,a,ja,ia,ao,jao,iao,nnzdmax,d,jd,id,
     &          dsub,jdsub,nsubmax,lindx,xlindx,nnzlmax,lnz,
     &          xlnz,iw(1,1),iw(1,2),iwmax,iwork,iw(1,3),iw(1,4),
     &          xsuper,iw(1,5),tmpmax,tmpvec,wwm(1,2),cachsz,
     &          level,x,s,u,c,y,b,wwn(1,1),wwn(1,2),wwn(1,3),
     &          wwn(1,4),nnzemax,e,je,ie,wwm(1,3),wwn(1,5),wwn(1,6),
     &          wwn(1,7),wwn(1,8),wwn(1,9),wwn(1,10),wwn(1,11),
     &          wwn(1,12),wwn(1,13),wwn(1,14),wwm(1,1),small,ierr,
     &          maxit,timewd)
      return
      end
      subroutine slpfn(n,m,nnza,a,ja,ia,ao,jao,iao,nnzdmax,d,jd,id,
     &                dsub,jdsub,nsubmax,lindx,xlindx,nnzlmax,lnz,
     &                xlnz,invp,perm,iwmax,iwork,colcnt,snode,xsuper,
     &                split,tmpmax,tmpvec,newrhs,cachsz,level,x,s,u,
     &                c,y,b,r,z,w,q,nnzemax,e,je,ie,dy,dx,ds,dz,dw,dxdz,
     &                dsdw,xi,xinv,sinv,ww1,ww2,small,ierr,maxit,timewd)
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c Sparse implentation of LMS's interior point method via 
c    Ng-Peyton's sparse Cholesky factorization for sparse 
c    symmetric positive definite
c INPUT:
c     n -- the number of row in the coefficient matrix A'
c     m -- the number of column in the coefficient matrix A'
c     nnza -- the number of non-zero elements in A'
c     a -- an nnza-vector of non-zero values of the design 
c          matrix (A') stored in csr format 
c     ja -- an nnza-vector of indices of the non-zero elements of
c           the coefficient matrix 
c     ia -- an (n+1)-vector of pointers to the begining of each
c           row in a and ja
c     ao -- an nnza-vector of work space for the transpose of
c           the design matrix stored in csr format or the 
c           design matrix stored in csc format
c     jao -- an nnza-vector of work space for the indices of the 
c            transpose of the design matrix 
c     iao -- an (n+1)-vector of pointers to the begining of each
c            column in ao and jao
c     nnzdmax -- upper bound of the non-zero elements in AA'
c     d -- an nnzdmax-vector of non-zero values of AQ^(-1)
c     jd -- an nnzdmax-vector of indices in d
c     id -- an (m+1)-vector of pointers to the begining of each
c           row in d and jd
c     dsub -- the values of e excluding the diagonal elements
c     jdsub -- the indices to dsub
c     nsubmax -- upper bound of the dimension of lindx
c     lindx -- an nsub-vector of interger which contains, in 
c           column major oder, the row subscripts of the nonzero
c           entries in L in a compressed storage format
c     xlindx -- an (m+1)-vector of integer of pointers for lindx
c     nnzlmax -- the upper bound of the non-zero entries in 
c                L stored in lnz, including the diagonal entries
c     lnz -- First contains the non-zero entries of d; later
c            contains the entries of the Cholesky factor
c     xlnz -- column pointer for L stored in lnz
c     invp -- an n-vector of integer of inverse permutation 
c             vector
c     perm -- an n-vector of integer of permutation vector
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
c     newrhs -- extra work vector for right-hand side and 
c               solution
c     cachsz -- size of the cache (in kilobytes) on the target 
c               machine
c     level -- level of loop unrolling while performing numerical
c              factorization
c     x -- an n-vector, the initial feasible solution in the primal
c              that corresponds to the design matrix A'
c     s -- an n-vector 
c     u -- an n-vector of upper bound for x
c     c -- an n-vector, usually the "negative" of
c          the pseudo response
c     y -- an m-vector, the initial dual solution 
c     b -- an n-vector, usualy the rhs of the equality constraint
c          X'a = (1-tau)X'e in the rq setting
c     r -- an n-vector of residuals
c     z -- an n-vector of the dual slack variable
c     w -- an n-vector 
c     q -- an n-vector of work array containing the diagonal 
c          elements of the Q^(-1) matrix
c     nnzemax -- upper bound of the non-zero elements in AA'
c     e -- an nnzdmax-vector containing the non-zero entries of
c          AQ^(-1)A' stored in csr format
c     je -- an nnzemax-vector of indices for e
c     ie -- an (m+1)-vector of pointers to the begining of each
c           row in e and je
c     dy -- work array
c     dx -- work array
c     ds -- work array
c     dz -- work array
c     dw -- work array
c     dxdz -- work array
c     dsdw -- work arry
c     xi -- work array
c     xinv -- work array
c     sinv -- work array
c     ww1 -- work array
c     ww2 -- work array
c     small -- convergence tolerance for inetrior algorithm
c     ierr -- error flag
c       1 -- insufficient work space in call to extract
c       2 -- nnzd > nnzdmax
c       3 -- insufficient storage in iwork when calling ordmmd;
c       4 -- insufficient storage in iwork when calling sfinit;
c       5 -- nnzl > nnzlmax when calling sfinit
c       6 -- nsub > nsubmax when calling sfinit
c       7 -- insufficient work space in iwork when calling symfct
c       8 -- inconsistancy in input when calling symfct
c       9 -- tmpsiz > tmpmax when calling bfinit; increase tmpmax
c       10 -- nonpositive diagonal encountered when calling 
c            blkfct, the matrix is not positive definite
c       11 -- insufficient work storage in tmpvec when calling 
c            blkfct
c       12 -- insufficient work storage in iwork when calling 
c            blkfct
c     maxit -- upper limit of the iteration; on return holds the
c              number of iterations
c     timew -- amount of time to execute this subroutine
c OUTPUT:
c     y -- an m-vector of primal solution
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
      integer nnza,m,n,nsuper,nnzdmax,nnzemax,iwmax,nnzd,
     &        nnzlmax,nsubmax,cachsz,level,tmpmax,ierr,maxit,it,
     &        ja(nnza),jao(nnza),jdsub(nnzemax+1),jd(nnzdmax),
     &        ia(n+1),iao(m+1),id(m+1),lindx(nsubmax),xlindx(m+1),
     &        invp(m),perm(m),xlnz(m+1),iwork(iwmax),
     &        colcnt(m),snode(m),xsuper(m+1),split(m),je(nnzemax),
     &        ie(m+1)
      double precision ddot,gap,zero,one,beta,small,deltap,deltad,mu,g,
     &                 a(nnza),ao(nnza),dsub(nnzemax+1),d(nnzdmax),
     &                 lnz(nnzlmax),c(n),b(m),newrhs(m),y(m),
     &                 tmpvec(tmpmax),r(n),z(n),w(n),x(n),s(n),
     &                 u(n),q(n),
     &                 e(nnzemax),dy(m),dx(n),ds(n),dz(n),dw(n),
     &                 dxdz(n),dsdw(n),xinv(n),sinv(n),xi(n),
     &                 ww1(n),ww2(m)
      double precision timewd(7)
      real gtimer,timbeg,timend
      external smxpy1,smxpy2,smxpy4,smxpy8
      external mmpy1,mmpy2,mmpy4,mmpy8
      parameter (beta=9.995d-1, one=1.0d0, zero=0.0d0)
      do i = 1,7
         timewd(i) = 0.0
      enddo
      it = 0
      nnzd = ie(m+1) - 1
      nnzdsub = nnzd - m
c
c  Compute the initial gap
c
       gap = ddot(n,z,1,x,1) + ddot(n,w,1,s,1)
c
c  Start iteration
c
   20 continue
      if(gap .lt. small .or. it .gt. maxit) goto 30
      it = it + 1
c
c  Create the diagonal matrix Q^(-1) stored in q as an n-vector
c  and update the residuals in r
c
      do i=1,n
         q(i) = one/(z(i)/x(i)+w(i)/s(i))
         r(i) = z(i) - w(i)
      enddo
c
c  Obtain AQ^(-1) and store in d,jd,id in csr format
c
      call amudia(m,1,ao,jao,iao,q,d,jd,id)
c
c  Obtain AQ^(-1)A' and store in e,je,ie in csr format
c
      call amub(m,m,1,d,jd,id,a,ja,ia,e,je,ie,nnzemax,iwork,ierr)
      if (ierr .ne. 0) then
         ierr = 2
         go to 100
      endif
c
c  Extract the non-diagonal structure of e,je,ie and store in dsub,jdsub
c
      call extract(e,je,ie,dsub,jdsub,m,nnzemax,nnzemax+1,ierr)
      if (ierr .ne. 0) then
         ierr = 1
         go to 100
      endif
c
c Compute b - Ax + AQ^(-1)r and store it in c in two steps
c  First: store Ax in ww2
      call amux(m,x,ww2,ao,jao,iao)
c
c  Second: save AQ^(-1)r in c temporarily
c
      call amux(m,r,c,d,jd,id)
      do i = 1,m
         c(i) = b(i) - ww2(i) + c(i)
      enddo
c
c  Compute dy = (AQ^(-1)A')^(-1)(b-Ax+AQ^(-1)r); result returned via dy
c
c Call chlfct to perform Cholesky's decomposition of e,je,ie
c
      call chlfct(m,xlindx,lindx,invp,perm,iwork,nnzdsub,jdsub,
     &            colcnt,nsuper,snode,xsuper,nnzlmax,nsubmax,xlnz,lnz,
     &            ie,je,e,cachsz,tmpmax,level,tmpvec,split,ierr,it,
     &            timewd)
      if (ierr .ne. 0) go to 100
c
c Call blkslv: Numerical solution for the new rhs stored in b
c
      do i = 1,m
         newrhs(i) = c(perm(i))
      enddo
      timbeg = gtimer()
      call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
      timend = gtimer()
      timewd(7) = timewd(7) + timend - timbeg
      do i = 1,m
         dy(i) = newrhs(invp(i))
      enddo
c
c  Compute dx = Q^(-1)(A'dy - r), ds = -dx, dz  and dw
c
      call amux(n,dy,dx,a,ja,ia)
      do i=1,n
         dx(i) = q(i) * (dx(i) - r(i))
         ds(i) = -dx(i)
         dz(i) = -z(i) * (one + dx(i) / x(i))
         dw(i) = -w(i) * (one + ds(i) / s(i))
      enddo
c
c Compute the maximum allowable step lengths
c
      call bound(x,dx,s,ds,z,dz,w,dw,n,beta,deltap,deltad)
      if (deltap * deltad .lt. one) then
c
c Update mu
c
         mu = ddot(n,z,1,x,1) + ddot(n,w,1,s,1)
         g = ddot(n,z,1,x,1) + deltap*ddot(n,z,1,dx,1)
     &       + deltad*ddot(n,dz,1,x,1) + deltad*deltap*ddot(n,dz,1,dx,1)
     &       + ddot(n,w,1,s,1) + deltap*ddot(n,w,1,ds,1)
     &       + deltad*ddot(n,dw,1,s,1) + deltad*deltap*ddot(n,dw,1,ds,1)
         mu = mu*((g/mu)**3)/(2.d0*dfloat(n))
c
c Compute dxdz and dsdw
c
         do i = 1,n
             dxdz(i) = dx(i)*dz(i)
             dsdw(i) = ds(i)*dw(i)
             xinv(i) = one/x(i)
             sinv(i) = one/s(i)
             xi(i) = xinv(i) * dxdz(i) - sinv(i) * dsdw(i)
     &               - mu * (xinv(i) - sinv(i))
             ww1(i) = q(i) * xi(i)
         enddo
c
c Compute AQ^(-1)(dxdz - dsdw - mu(X^(-1) - S^(-1))) and
c store it in ww2 temporarily
c
         call amux(m,ww1,ww2,ao,jao,iao)
         do i = 1,m
            c(i) = c(i) + ww2(i)
         enddo
c
c
c Compute dy and return the result in dy
c
c Call blkslv: Numerical solution for the new rhs stored in b
c
      do i = 1,m
         newrhs(i) = c(perm(i))
      enddo
      timbeg = gtimer()
      call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
      timend = gtimer()
      timewd(7) = timewd(7) + timend - timbeg
      do i = 1,m
         dy(i) = newrhs(invp(i))
      enddo
c
c  Compute dx = Q^(-1)(A'dy - r + mu(X^(-1) - S^(-1)) -dxdz + dsdw),
c  ds = -dx, dz  and dw
c
         call amux(n,dy,dx,a,ja,ia)
         do i=1,n
            dx(i) = q(i) * (dx(i) - xi(i) - r(i))
            ds(i) = -dx(i)
            dz(i) = -z(i) + xinv(i)*(mu - z(i)*dx(i) - dxdz(i))
            dw(i) = -w(i) + sinv(i)*(mu - w(i)*ds(i) - dsdw(i))
         enddo
c
c Compute the maximum allowable step lengths
c
         call bound(x,dx,s,ds,z,dz,w,dw,n,beta,deltap,deltad)
      endif
c
c Take the step
c
      call daxpy(n,deltap,dx,1,x,1)
      call daxpy(n,deltap,ds,1,s,1)
      call daxpy(n,deltad,dw,1,w,1)
      call daxpy(n,deltad,dz,1,z,1)
      call daxpy(m,deltad,dy,1,y,1)
      gap = ddot(n,z,1,x,1) + ddot(n,w,1,s,1)
      goto 20
   30 continue
  100 continue
      maxit = it
      return
      end
