C  Comments, bug reports, etc. should be sent to:  portnoy@stat.uiuc.edu
C
      SUBROUTINE CRQ(M,N,MPLUS,N2,X,Y,C,TZERO,MAXW,STEP,IFT,H,XH,
     * WA,WB,WC,WD,WE,WF,IA,NSOL,SOL,LSOL,ICEN,TCEN,LCEN)
C
C     M = number of Observations
C     N = Number of Parameters
C     MPLUS = row dim of X, WA, and WC .GE. M+1
C     N2 = N+2
C     X is the X matrix (MPLUS BY N) (includes censored obs)
C     Y is the y vector (M) (includes censored obs)
C     C is the censoring integer(M) vector, 1 if obs is censored,
C                                           0 otherwise
C     TZERO is the initial tau value
C     MAXW ( .LE. MPLUS ): Maximal number of calls to weighted RQ
C           for possible degeneracies. If exceeded (IFT = 7), dither.
C           If  MAXW .LE. 0, don't pivot -- use only weighted rq
C     STEP: Step size for weighted rq at degeneracy
C     IFT exit code:
C        0:  ok
C        1:  M < 0  or  N < 0  OR  M < N
C	 2:  MPLUS .LT. M + 1  or  N2 .NE. N+2
C        3:  MAXW > MPLUS
C        4:  less than N noncensored obs above the tau=0 soln
C        5:  possible degeneracy, rq called at tau + step, see IA
C        6:  LSOL becomes greater than NSOL
C        7:  MAXW exceeded
C        8:  weighted rq tries to include infinite basis element
C     H  is an integer work vector (N) ( = basis indices )
C     XH is a double precision work array (N by N)
C     WA is a double precision work array (MPLUS by N)
C     WB is a double precision work vector (MPLUS)
C     WC is a double precision work array (MPLUS by N+2)
C     WD is a double precision work vector (MPLUS)
C     WE is a double precision work vector (MPLUS)
C     WF is a double precision work vector (N)
C     IA is an integer work vector (MPLUS) 
C     NSOL is an (estimated) row dimension of the solution array
C          if NSOL < M, solutions returned only at tau = i/(NSOL-1)
C          if all solutions are desired, NSOL = 3*M is a good choice
C  on output:
C     SOL is the coefficient solution array (N+2 by NSOL)
C         first and last rows give tau bkpts and quantile
C         rows 2:(N+1) give the beta coefficients
C     LSOL is the actual dimension of the solution arrays
C          LSOL = NSOL  if NSOL < M, else LSOL .LE. NSOL
C     ICEN (M) indicates status of censored observations
C          = 0 if uncensored (or uncrossed censored while pivoting)
C          = 1 if crossed censored
C          = 2 if deleted as below tau = 0 solution
C          = 3 if above last (maximal tau) solution
C     TCEN (M) are the corresponding tau value where censored obs is
C          firstcrossed:  weight = (tau - TCEN(I))/(1 - TCEN(I))
C          TCEN = 1 if censored obs is never crossed (ICEN = 3)
C          TCEN = 0 if censored obs is deleted (ICEN = 2)
C     LCEN = number of censored observations
C     WD = list of first MPLUS tau values at which degeneracy
C          may have occurred (and weighted RQ was called)
C     H(1) = number of calls to weighted RQ (apparent degeneracy)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER H(N),ICEN(M),C(M),IA(MPLUS)
      LOGICAL APC, DEG
      DOUBLE PRECISION ONE
      DIMENSION X(M,N),Y(M),SOL(N2,NSOL),WC(MPLUS,N2),XH(N,N)
      DIMENSION WA(MPLUS,N),WB(MPLUS),WD(MPLUS),WE(MPLUS)
      DIMENSION F(2),WF(N),TCEN(M)
      DATA BIG/1.D17/
      DATA ZERO/0.00D0/
      DATA HALF/0.50D0/
      DATA ONE/1.00D0/
      DATA TWO/2.00D0/
C
C  CHECK DIMENSION PARAMETERS
C
      DEG = .FALSE.
      IFT=0
      N1 = NSOL-1
      LCEN = 0
      DO 2 I = 1,M
      IF(C(I) .EQ. 1) LCEN = LCEN + 1
   2  CONTINUE
      IF(MPLUS .LT. M+1) IFT = 2
      IF(N2 .NE. N+2) IFT = 2
      IF(M .LE. ZERO .OR. N. LE. ZERO .OR. M .LT. N) IFT = 1
      IF(MAXW .GT. MPLUS) IFT = 3
      IF(IFT .NE. 0) RETURN
C
C  INITIALIZATION
C
      TOLER = 1.D-11
      TOL1 = 10.0D0*TOLER
      IF(TZERO .LE. ZERO  .OR.  TZERO .GT. STEP) TZERO = STEP
      DO 5 I=1,M
        ICEN(I) = 0
        TCEN(I) = ONE
        WB(I) = Y(I)
        DO 5 J=1,N
   5    WA(I,J) = X(I,J)
      MA = M
C
C  GET TAU = O SOLUTION
C  DELETE CENSORED OBS BELOW TAU = 0 SOLN
C
  10  IF(MA .LE. N) THEN
        IFT = 4
        RETURN
      ENDIF
      CALL RQ1(MA,N,MPLUS,N2,WA,WB,TZERO,TOLER,IFT1,
     * WF,WE,IA,WC,WD)
      K = 0
      KL = 0
C
C  CHECK FOR CENSORED OBS BELOW TAU = 0 SOLN
C
      L = 0
      DO 20 I=1,M
      IF(ICEN(I) .EQ. 2) GO TO 20
      L = L + 1
      IF( C(I) .EQ. 1 .AND. WE(L) .LE. TOL1 ) THEN
        K = K + 1
        ICEN(I) = 2
        GO TO 20
      ENDIF
      IF( DABS(WE(L)) .GE. TOL1) GO TO 20
      KL = KL + 1
      IF(KL .LE. N) H(KL) = I
  20  CONTINUE
      IF( K .EQ. 0) GO TO 40
      MA = MA - K
      L = 0
C
C  DELETE CENSORED OBS BELOW TAU = 0 SOLN, AND RETRY
C 
      DO 30 I=1,M
        IF(ICEN(I) .EQ. 2) GO TO 30
        L = L + 1
        WB(L) = Y(I)
        DO 25 J=1,N
  25      WA(L,J) = X(I,J)
  30    CONTINUE
      GO TO 10
  40  CONTINUE
C
C  SAVE INITIAL SOLUTION, AND INITIALIZE PIVOTING LOOP
C
      SOL(1,1) = ZERO
      DO 50 J = 1,N
  50  SOL(J+1,1) = WF(J)
      NWRQ = 0
      APC = .FALSE.
      LSOL = 2
      TAU = TZERO
C
C  COMPUTE NEXT TAU
C  FIRST CHECK FOR DEGENERACY AT TZERO
C
      IF (KL .GT. N)  GO TO 500
C
C  GET X(H,)
C
      DO 70 K = 1,N
        I= H(K)
        DO 60 J = 1,N
  60    XH(K,J) = X(I,J)
  70  CONTINUE
C
C  GET X(H,) INVERSE
C
      CALL DGECO(XH,N,N,IA,V,WC)
      IF(V. LT. TOLER) GO TO 500
      JOB = 01
      CALL DGEDI(XH,N,N,IA,F,WC,JOB)
C
C  GET RESIDUALS
C
      DO 90 I=1,M
      S = Y(I)
      DO 80 J = 1,N
  80  S = S - WF(J)*X(I,J)
  90  WE(I) = S
C
C  BEGIN PIVOTING; WD = GRAD UPPER BD, IA = SIGN(GRAD DENOM)
C
 200  CONTINUE
      CALL GRAD(X,Y,M,N,H,ICEN,TCEN,XH,WE,TOLER,IA,WC,WD)
      KL = 0
      KM = 1
      S = WD(1)
      IF(N .EQ. 1) GO TO 230
      DO 210 J=2,N
 210  S = DMIN1(WD(J),S)
      DO 220 J=1,N
      IF(WD(J) .GE. S + TOLER) GO TO 220
        KL = KL + 1
        KM = J
 220  CONTINUE
 230  CONTINUE
C     
C  IF ALL POS RESIDS CENSORED, RETURN; ELSE CHECK FOR DEGEN.
C     
      IF(APC) THEN
        SOL(1,LSOL) = DMAX1(S, TAU)
        DO 240 J=1,N
 240    SOL(J+1,LSOL) = WF(J)
        GO TO 600
      ENDIF
C
C  IF DEGENERACY, CALL WEIGHTED RQ
C
      IF (KL .GT. 1) GO TO 500
C     
C  CHECK FOR INFEASIBILITY CAUSED BY REWEIGHTING
C 
      IF(S .LT. TAU + TOL1) THEN
        LSOL = LSOL - 1
        TAU = TAU - .8*STEP
        GO TO 500
      ENDIF
      TAU = S
C
C  GET NEW OBSERVATION TO ENTER BASIS
C  FIRST DEFINE WD = BASIS INDICATOR = 1 IF I IN H(J)
C
      DO 250 I=1,M
 250  WD(I) = ZERO
      DO 260 J=1,N
 260  WD(H(J)) = ONE
      KN = 0
      D = BIG 
      KIN = 0
      DO 300 I=1,M
      IF(ICEN(I) .EQ. 2) GO TO 300
      S = ZERO
      DO 270 J=1,N
 270  S = S + X(I,J)*XH(J,KM) 
      S = IA(KM)*S
      IF(DABS(S).LT.TOL1 .OR. (C(I).EQ.1 .AND. ICEN(I).NE.1)
     *    .OR. WD(I).EQ.ONE) GO TO 300
      S = WE(I)/S
      IF(S .LT. TOL1) GO TO 300
      IF(S .LE. D - TOL1) THEN
        D = S
        KIN = I
        KN = 0
      ELSE
        IF(S .LE. D + TOL1) KN = 1
      ENDIF
 300  CONTINUE
      IF(KN .EQ. 1) GO TO 500
      H(KM) = KIN
C
C  UPDATE SOLUTION
C  GET NEW X(H,)
C
        DO 310 K = 1,N
        I = H(K)
        DO 310 J = 1,N
 310    XH(K,J) = X(I,J)
C
C  GET X(H,) INVERSE
C
 315  CALL DGECO(XH,N,N,IA,V,WA)
      IF(V. LT. TOLER) GO TO 500
      JOB = 01
      CALL DGEDI(XH,N,N,IA,F,WA,JOB)
C
C  GET BETA-HAT, RESIDS
C
      DO 340 K = 1,N
      S = ZERO
      DO 330 J = 1,N
 330  S = S + XH(K,J)*Y(H(J))
 340  WF(K) = S
 345  DO 360 I=1,M
      S = Y(I)
      DO 350 J = 1,N
 350  S = S - WF(J)*X(I,J)
 360  WE(I) = S
C
C  SAVE SOLUTION
C
 365  DO 380 I = (LSOL-1),N1
      IF(NSOL.GE.M .OR. TAU .GT. DFLOAT(I)/DFLOAT(N1) - 10*TOL1) THEN
        SOL(1,LSOL) = TAU
        DO 370 J = 1,N
 370    SOL(J+1,LSOL) = WF(J)
        LSOL = LSOL + 1
        GO TO 390
      ENDIF
 380  CONTINUE
 390  CONTINUE
C
C  CHECK FOR DIM(SOL) EXCEEDED
C
      IF(LSOL .GT. NSOL) THEN
        IFT = 6
        GO TO 600
      ENDIF
      IF(APC) THEN
        LSOL = LSOL - 1
        GO TO 600
      ENDIF
C
C  SET WT FOR CROSSED CENSORED OBSERVATIONS
C  APC = .TRUE.  IF ALL POS RESID ARE UNCROSSED CENSORED
C
      APC = .TRUE.
C  if the following are left uncommented: still get zero column
C      TAUW = TAU - STEP/TWO
C      IF(MAXW.GT.ZERO)THEN
C        TAUW=TAU
C      ENDIF
      DO 400 I=1,M
      IF(ICEN(I) .EQ. 2) GO TO 400
      IF(WE(I).GE.TOL1 .AND. (C(I).EQ.0 .OR.  ICEN(I).EQ.1))
     *   APC = .FALSE.
      IF(WE(I).LE.-TOL1 .AND. C(I).EQ.1 .AND. ICEN(I).EQ.0) THEN
C  
C at this point,  Y(I) is a crossed censored obs  (C_i)
C for the grid method  TAUW = TAU - STEP*(B2-Y(I))/(B2-B1)
C   where   B1 = x_i' beta_j   and   B2 = x_i' beta_(j+1)
C otherwise (for pivot), TAUW = TAU
C
        IF(MAXW.LT.0) THEN
         B1 = ZERO
         B2 = ZERO
         DO 396 J=1,N
          B1 = B1 + X(I,J) * SOL(J+1,LSOL - 2)
396       B2 = B2 + X(I,J) * SOL(J+1,LSOL - 1)
         A1 = (B2 - Y(I))/(B2-B1)
         TAUW = TAU - A1*STEP
        ELSE
         TAUW = TAU
        ENDIF
        ICEN(I) = 1
        TCEN(I) = TAUW
      ENDIF
 400  CONTINUE
      IF(APC) THEN
        IF(DEG) THEN
          IFT = 5
          LSOL = LSOL - 1
          GO TO 600
        ENDIF
        GO TO 200
      ENDIF
C
C  RETURN IF TAU > 1
C
      IF(TAU .GE. ONE - 10.0D0 * TOL1) GO TO 600
      IF(DEG) GO TO 500
      IF(MAXW .GT. 0) GO TO 200
C
C  POSSIBLE DEGENERACY -- USE WEIGHTED RQ1 TO RESOLVE
C
 500  NWRQ = NWRQ + 1
      IF(MAXW .GT. 0  .AND.  NWRQ .GT. MAXW) THEN
        IFT = 7
        LSOL = LSOL - 1
        GO TO 600
      ENDIF
      IF(NWRQ .LE. NSOL) SOL(N+2,NWRQ) = TAU
      TAU = TAU + STEP
      IF(TAU .GE. ONE) TAU = ONE - TOL1
      L = 0
      L1 = 0
      DO 506 J=1,N
 506  WF(J) = ZERO
      DO 530 I = 1,M
      IF(ICEN(I) .EQ. 0) THEN
        IF (C(I) .EQ. 0) THEN
          L1 = L1 + 1
          WB(L1) = Y(I)
          DO 510 J=1,N
 510      WA(L1,J) = X(I,J)
        ELSE
          DO 514 J = 1,N
 514      WF(J) = WF(J) + X(I,J)
        ENDIF
      ENDIF
      IF(ICEN(I) .EQ. 1) THEN
        L1 = L1 + 1
        L = L+1
        W = (TAU - TCEN(I))/(ONE - TCEN(I))
        WB(L1) = W * Y(I)
        DO 520 J = 1,N
          WA(L1,J) = W * X(I,J)
 520      WF(J) = WF(J) + (ONE - W) * X(I,J)
      ENDIF
 530  CONTINUE
      MAL = L1+1
      DO 534 J=1,N
 534  WA(MAL,J) = WF(J)
      WB(MAL) = BIG
      CALL RQ1(MAL,N,MPLUS,N2,WA,WB,TAU,TOLER,IFT1,WF,WE,IA,WC,WD)
      DEG = .FALSE.
      IF(DABS(WE(MAL)) .LE. .0001) THEN
        IFT = 8
        LSOL = LSOL - 1
        GO TO 600
      ENDIF 
      L = 0
      K = 0
      DO 550 I=1,M
        IF(ICEN(I) .EQ. 2) GO TO 550
        L = L+1
        IF(DABS(WE(L)) .GE. TOL1) GO TO 550
          K = K+1
          IF(K .LE. N) THEN
            H(K) = I
            DO 540 J=1,N
 540        XH(K,J) = X(I,J)
          ENDIF
 550  CONTINUE
      IF(K .LT. N) THEN
        IFT = 8
        LSOL = LSOL - 1
        GO TO 600
      ENDIF
      IF(K .GT. N) DEG = .TRUE.
      GO TO 345
C
C  DEFINE OUTPUT AND RETURN
C
 600  H(1) = NWRQ
      L = MIN0(NWRQ,MPLUS)
      DO 610 I=1,L
 610  WD(I) = SOL(N+2,I)
      DO 620 I=1,M
      IF(ICEN(I) .EQ. 2) GO TO 620
      IF(C(I) .EQ.1 .AND. TCEN(I) .EQ. ONE) ICEN(I) = 3
 620  CONTINUE
      V = DFLOAT(MA)
      DO 650 J = 1,N
      S = ZERO
        DO 640 I = 1,M
          IF(ICEN(I) .EQ. 2) GO TO 640
          S = S + X(I,J)
 640    CONTINUE
 650  WF(J) = S/V
      DO 670 I = 1,LSOL
      S = ZERO
        DO 660 J = 1,N
 660    S = S + WF(J) * SOL(J+1,I)
 670  SOL(N+2,I) = S
      RETURN
      END
      SUBROUTINE GRAD(X,Y,M,N,H,ICEN,TCEN,XH,
     * R,TOL,IFLAG,WA,GUP)
C
C  X matrix (M BY N)
C  Y vector (M)
C  M = Number of Observations
C  N = Number of Parameters
C  H = basis, integer(N) vector 
C  ICEN (integer(M)) = 2 for deleted obs
C                    = 1 for crossed non-censored obs
C                    = 0 otherwise
C  TCEN(M) are the corresponding tau values where a censored obs is
C    first crossed (ICEN=1): weight = (tau - TCEN(I))/(1 - TCEN(I))
C    TCEN = 1 if censored obs is never crossed
C  XH = (N by N) X(H,) inverse matrix
C  R(M) = residuals
C  TOL = zero tolerance (1.D-10 by default)
C  IFLAG (M)  work vector
C  WA (M by N) work array
C  returns:
C    GUP(N)  = upper bounds on tau
C    IFLAG(1:N) = sign(grad denom)
C               = +1 if bound from lower grad bound
C               = -1 if bound from upper grad bound
C    
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER H(N),ICEN(M),IFLAG(*)
      DOUBLE PRECISION ONE
      DIMENSION X(M,N),Y(M),GUP(N),XH(N,N)
      DIMENSION R(M),TCEN(M),WA(M,N)
      DATA  ZERO/0.00D0/, ONE/1.0d0/
C
C GET X'(XH-INVERSE)
C
      DO 80 I = 1,M
      IF(ICEN(I) .EQ. 2) GO TO 80
      DO 70 J = 1,N
      A = ZERO
        DO 60 K = 1,N
  60    A = A + X(I,K)*XH(K,J)
  70  WA(I,J) = A
  80  CONTINUE
  
C
C  GET GRADIENT BOUNDS
C  FIRST SET IFLAG = 1 FOR BASIS INDICES (TEMPORARY)
C
      T = ZERO
      DO 90 I = 1,M
  90  IFLAG(I) = 0
      DO 95 J = 1,N
  95  IFLAG(H(J)) = 1
      DO 120 J=1,N
      A = ZERO
      B = ZERO
      C = ZERO
      D = ZERO
      DO 100 I = 1,M
        IF(ICEN(I) .EQ. 2) GO TO 100
        IF(ICEN(I) .EQ. 0) THEN
          IF(R(I) .GT. TOL)    A = A + WA(I,J)
          IF(R(I) .LT. - TOL)  B = B + WA(I,J)
          GO TO 100
        ENDIF
        IF(IFLAG(I).EQ.1) GO TO 100
C
C  IFLAG = 0: NOT BASIS  AND  ICEN = 1: CROSSED CENSORED
C
        IF(R(I) .LT. -TOL) THEN
           T = TCEN(I)/(ONE - TCEN(I))
           C = C - T*WA(I,J)
           GO TO 100
        ENDIF
        IF(R(I) .GT. TOL) THEN          
           D = D - WA(I,J)
        ENDIF
 100  CONTINUE
      D = D - C
      L = ICEN(H(J))
      IF(L .NE. 0) T = TCEN(H(J))/( ONE - TCEN(H(J)) )
      S = DFLOAT(L)*(T + ONE) - ONE
      E1 = A + B - D - S
      E2 = A + B - D + ONE
      IF(E1 .GT. ZERO) THEN
        GUP(J) = (B + C - S)/E1
        IFLAG(J+M) = 1
        GO TO 120
      ENDIF
      IF(E2 .LT. ZERO) THEN
        GUP(J) = (B + C)/E2
        IFLAG(J+M) = -1
        GO TO 120
      ENDIF
      GUP(J) = - ONE
 120  CONTINUE
      DO 130 J = 1,N
 130  IFLAG(J) = IFLAG(J+M)
      RETURN
      END
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,*),det(2),work(*)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(*)
      double precision a(lda,*),z(*)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
c      use numerical_libraries
	
      integer lda,n,ipvt(*),info
      double precision a(lda,*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
c           integer idamax,j,k,kp1,l,nm1
	 integer j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine  dscal(n,da,dx,incx)
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end



      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
