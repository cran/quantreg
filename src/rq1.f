C	This is a very primative version of the barrodale-roberts rq
C	algorithm intended only to be used for resampling procedures.
C       bug reports should be sent to roger@ysidro.econ.uiuc.edu or to
C
C				Roger Koenker
C				Department of Economics
C				University of Illinois
C				Champaign, Illinois, 61820
C
C
      SUBROUTINE RQ1(M,N,M5,N2,A,B,TAU,TOLER,IFT,X,E,S,WA,WB,
     *NSOL,SOL,DSOL,LSOL)
C
C     Number of Observations
C     Number of Parameters
C     M+5
C     N+2
C     A is the X matrix
C     B is the y vector
C     tau, the desired quantile
C     toler, smallest detectable |x-y|/x machine precision to the 2/3
C     ift exit code:
C		0-ok
C		else dimensions inconsistent
C     X the parameter estimate betahat
C     E is the residual vector
C     S is an integer work array
C     WA is a work array
C     WB is another work array
C     NSOL is an estimated (row) dimension of the solution array say 3*M
C     SOL is the primal solution array
C     DSOL is the dual solution arry
C     LSOL is the actual dimension of the solution arrays
C     Utilization:  If you just want a solution at a single quantile you
C     needn't bother with SOL, NSOL, etc, if you want all the solutions
C     then set TAU to something <0 and SOL and DSOL will return all the
C     estimated quantile solutions.
C     The algorithm  is a slightly modified version of Algorithm AS 229 
C     described in Koenker and D'Orey, "Computing Regression Quantiles,
C     Applied Statistics, pp. 383-393. 
C
      INTEGER I,J,K,KL,KOUNT,KR,L,LSOL,M,M1,M2,M3,M4,M5,IFT
      INTEGER N,N1,N2,NSOL,OUT,S(M)
      LOGICAL STAGE,TEST,INIT,IEND
      DOUBLE PRECISION A1,AUX,B1,BIG,D,DIF,PIVOT,SMAX,TAU,T0,T1,TNT
      DOUBLE PRECISION MIN,MAX,TOLER,ZERO,HALF,ONE,TWO
      DOUBLE PRECISION B(M),SOL(N2,NSOL),A(M,N),X(N),WA(M5,N2),WB(M)
      DOUBLE PRECISION SUM,E(M),DSOL(M,NSOL)
      PARAMETER( BIG = 1.D37)
      PARAMETER( ZERO = 0.D00)
      PARAMETER( HALF = 0.5D0)
      PARAMETER( ONE = 1.D0)
      PARAMETER( TWO = 2.D0)
C
C  CHECK DIMENSION PARAMETERS
C
      IFT=0
      IF(M5 .NE. M+5)IFT = 3
      IF(N2 .NE. N+2)IFT = 4
      IF(M.LE.ZERO.OR.N.LE.ZERO)IFT = 5
      IF(IFT .GT. TWO)RETURN
C
C  INITIALIZATION
C
      M1 = M+1
      N1 = N+1
      M2 = M+2
      M3 = M+3
      M4 = M+4
      DO 2 I=1,M
      WB(I)=B(I)
      DO 1 J=1,N
      WA(I,J)=A(I,J)
  1   CONTINUE
  2   CONTINUE
      WA(M2,N1)=ZERO
      DIF = ZERO
      IEND = .TRUE.
      IF(TAU .GE. ZERO .AND. TAU .LE. ONE)GOTO 3
      T0 = ONE/FLOAT(M)-TOLER
      T1 = ONE - T0
      TAU = T0
      IEND = .FALSE.
  3   CONTINUE
      INIT = .FALSE.
      LSOL = 1
      KOUNT = 0
      DO 9 K=1,N
      WA(M5,K) = ZERO
      DO 8 I=1,M
      WA(M5,K) = WA(M5,K) + WA(I,K)
  8   CONTINUE
      WA(M5,K) = WA(M5,K)/FLOAT(M)
  9   CONTINUE
      DO 10 J=1,N
      WA(M4,J) = J
      X(J) = ZERO
 10   CONTINUE
      DO 40 I=1,M
      WA(I,N2) = N+I
      WA(I,N1) = WB(I)
      IF(WB(I).GE.ZERO)GOTO 30
      DO 20 J=1,N2
      WA(I,J) = -WA(I,J)
 20   CONTINUE
 30   E(I) = ZERO
 40   CONTINUE
      DO 42 J=1,N
      WA(M2,J) = ZERO
      WA(M3,J) = ZERO
      DO 41 I=1,M
      AUX = SIGN(ONE,WA(M4,J)) * WA(I,J)
      WA(M2,J) = WA(M2,J) + AUX * (ONE - SIGN(ONE,WA(I,N2)))
      WA(M3,J) = WA(M3,J) + AUX * SIGN(ONE,WA(I,N2))
 41   CONTINUE
      WA(M3,J) = TWO * WA(M3,J)
 42   CONTINUE
      GOTO 48
 43   CONTINUE
      LSOL = LSOL + 1
      DO 44 I=1,M
      S(I) = ZERO
 44   CONTINUE
      DO 45 J=1,N
      X(J) = ZERO
 45   CONTINUE
C
C  COMPUTE NEXT TAU
C
      SMAX = TWO
      DO 47 J=1,N
      B1 =  WA(M3,J)
      A1 = (-TWO - WA(M2,J))/B1
      B1 = -WA(M2,J)/B1
      IF(A1 .LT. TAU)GOTO 46
      IF(A1 .GE. SMAX) GOTO 46
      SMAX = A1
      DIF = (B1 - A1 )/TWO
 46   IF(B1 .LE. TAU) GOTO 47
      IF(B1 .GE. SMAX)GOTO 47
      SMAX = B1
      DIF = (B1 - A1)/TWO
 47   CONTINUE
      TNT = SMAX + TOLER * (ONE + ABS(DIF)) 
      IF(TNT .GE. T1 + TOLER)IEND = .TRUE.
      TAU = TNT
      IF(IEND)TAU = T1
 48   CONTINUE
C
C COMPUTE NEW MARGINAL COSTS
C
      DO 49 J=1,N
      WA(M1,J) = WA(M2,J) + WA(M3,J) * TAU
 49   CONTINUE
      IF(INIT) GOTO 265
C
C STAGE 1
C
C DETERMINE THE VECTOR TO ENTER THE BASIS
C
      STAGE=.TRUE.
      KR=1
      KL=1
 70   MAX=-ONE
      DO 80 J=KR,N
      IF(ABS(WA(M4,J)).GT.N)GOTO 80
      D=ABS(WA(M1,J))
      IF(D.LE.MAX)GOTO 80
      MAX=D
      IN=J
 80   CONTINUE
      IF(WA(M1,IN).GE.ZERO)GOTO 100
      DO 90 I=1,M4
      WA(I,IN)=-WA(I,IN)
 90   CONTINUE
C
C DETERMINE THE VECTOR TO LEAVE THE BASIS
C
 100  K=0
      DO 110 I=KL,M
      D=WA(I,IN)
      IF(D.LE.TOLER)GOTO 110
      K=K+1
      WB(K)=WA(I,N1)/D
      S(K)=I
      TEST=.TRUE.
 110  CONTINUE
 120  IF(K.GT.0)GOTO 130
      TEST=.FALSE.
      GOTO 150
 130  MIN=BIG
      DO 140 I=1,K
      IF(WB(I).GE.MIN)GOTO 140
      J=I
      MIN=WB(I)
      OUT=S(I)
 140  CONTINUE
      WB(J)=WB(K)
      S(J)=S(K)
      K=K-1
C
C CHECK FOR LINEAR DEPENDENCE IN STAGE 1
C
 150  IF(TEST.OR..NOT.STAGE)GOTO 170
      DO 160 I=1,M4
      D=WA(I,KR)
      WA(I,KR)=WA(I,IN)
      WA(I,IN)=D
 160  CONTINUE
      KR=KR+1
      GOTO 260
 170  IF(TEST)GOTO 180
      WA(M2,N1)=TWO
      GOTO 390
 180  PIVOT=WA(OUT,IN)
      IF(WA(M1,IN)-PIVOT-PIVOT.LE.TOLER)GOTO 200
      DO 190 J=KR,N1
      D=WA(OUT,J)
      WA(M1,J)=WA(M1,J)-D-D
      WA(M2,J)=WA(M2,J)-D-D
      WA(OUT,J)=-D
 190  CONTINUE
      WA(OUT,N2)=-WA(OUT,N2)
      GOTO 120
C
C PIVOT ON WA(OUT,IN)
C
 200  DO 210 J=KR,N1
      IF(J.EQ.IN)GOTO 210
      WA(OUT,J)=WA(OUT,J)/PIVOT
 210  CONTINUE
      DO 230 I=1,M3
      IF(I.EQ.OUT)GOTO 230
      D=WA(I,IN)
      DO 220 J=KR,N1
      IF(J.EQ.IN)GOTO 220
      WA(I,J)=WA(I,J)-D*WA(OUT,J)
 220  CONTINUE
 230  CONTINUE
      DO 240 I=1,M3
 
      IF(I.EQ.OUT)GOTO 240
      WA(I,IN)=-WA(I,IN)/PIVOT
 240  CONTINUE
      WA(OUT,IN)=ONE/PIVOT
      D=WA(OUT,N2)
      WA(OUT,N2)=WA(M4,IN)
      WA(M4,IN)=D
      KOUNT=KOUNT+1
      IF(.NOT.STAGE)GOTO 270
C
C INTERCHANGE ROWS IN STAGE 1
C
      KL=KL+1
      DO 250 J=KR,N2
      D=WA(OUT,J)
      WA(OUT,J)=WA(KOUNT,J)
      WA(KOUNT,J)=D
 250  CONTINUE
 260  IF(KOUNT+KR.NE.N1)GOTO 70
C
C STAGE 2
C
 265  STAGE=.FALSE.
C
C DETERMINE THE VECTOR TO ENTER THE BASIS
C
 270  MAX=-BIG
      DO 290 J=KR,N
      D=WA(M1,J)
      IF(D.GE.ZERO)GOTO 280
      IF(D.GT.(-TWO))GOTO 290
      D=-D-TWO
 280  IF(D.LE.MAX)GOTO 290
      MAX=D
      IN=J
 290  CONTINUE
      IF(MAX.LE.TOLER)GOTO 310
      IF(WA(M1,IN).GT.ZERO)GOTO 100
      DO 300 I=1,M4
      WA(I,IN)=-WA(I,IN)
 300  CONTINUE
      WA(M1,IN)=WA(M1,IN)-TWO
      WA(M2,IN)=WA(M2,IN)-TWO
      GOTO 100
C
C COMPUTE QUANTILES
C
 310  CONTINUE
      DO 320 I=1,KL-1
      K=WA(I,N2)*SIGN(ONE,WA(I,N2))
      X(K) = WA(I,N1) * SIGN(ONE,WA(I,N2))
 320  CONTINUE
      SUM=ZERO
      DO 330 I=1,N
      SUM=SUM + X(I) * WA(M5,I)
      SOL(I+2,LSOL)=X(I)
      KD=ABS(WA(M4,I))-N
      DSOL(KD,LSOL)=ONE+WA(M1,I)/TWO
      IF(WA(M4,I).LT.ZERO)DSOL(KD,LSOL)=ONE-DSOL(KD,LSOL)
 330  CONTINUE
      DO 335 I=KL,M
      KD=ABS(WA(I,N2))-N
      IF(WA(I,N2).LT.ZERO)DSOL(KD,LSOL)=ZERO
      IF(WA(I,N2).GT.ZERO)DSOL(KD,LSOL)=ONE
 335  CONTINUE
      SOL(1,LSOL) = SMAX
      SOL(2,LSOL) = SUM
      IF(IEND)GOTO 340
      INIT = .TRUE.
      GOTO 43
 340  CONTINUE
      IF(LSOL.LE.2)GOTO 355
      SOL(1,1)=ZERO
      SOL(1,LSOL)=ONE
      DO 350 I=1,M
      DSOL(I,1)=ONE
      DSOL(I,LSOL)=ZERO
 350  CONTINUE
 355  CONTINUE
      L =KL-1
      DO 370 I=1,L
      IF(WA(I,N1).GE.ZERO)GOTO 370
      DO 360 J=KR,N2
      WA(I,J)=-WA(I,J)
 360  CONTINUE
 370  CONTINUE
      WA(M2,N1)=ZERO
      IF(KR.NE.1)GOTO 390
      DO 380 J=1,N
      D=ABS(WA(M1,J))
      IF(D.LE.TOLER.OR.TWO-D.LE.TOLER)GOTO 390
 380  CONTINUE
      WA(M2,N1)=ONE
 390  SUM = ZERO
      DO 400 I=KL,M
      K = WA(I,N2) * SIGN(ONE,WA(I,N2))
      D = WA(I,N1) * SIGN(ONE,WA(I,N2))
      SUM = SUM + D * SIGN(ONE,D) * (HALF + SIGN(ONE,D)*(TAU-HALF))
      K=K-N
      E(K)=D
 400  CONTINUE
      WA(M2,N2)=KOUNT
      WA(M1,N2)=N1-KR
      WA(M1,N1) = SUM
      RETURN
      END
