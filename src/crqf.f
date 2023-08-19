C  Comments, bug reports, etc. should be sent to:  portnoy@stat.uiuc.edu
C
      SUBROUTINE CRQF(M,N,MPLUS,N2,X,Y,C,TZERO,MAXW,STEP,IFT,H,XH,
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
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      CALL GRAD(X,M,N,H,ICEN,TCEN,XH,WE,TOLER,IA,WC,WD)
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
      IF(NSOL.GE.M .OR. TAU .GT. DBLE(I)/DBLE(N1) - 10*TOL1) THEN
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
      V = DBLE(MA)
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
      SUBROUTINE GRAD(X,M,N,H,ICEN,TCEN,XH,
     * R,TOL,IFLAG,WA,GUP)
C
C  X matrix (M BY N)
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
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER H(N),ICEN(M),IFLAG(*)
      DOUBLE PRECISION ONE
      DIMENSION X(M,N),GUP(N),XH(N,N)
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
      S = DBLE(L)*(T + ONE) - ONE
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
