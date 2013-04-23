      SUBROUTINE SRTPAI ( A, SA, P, SP, N )
C
C     SRTPAI SETS P(1) = 1, P(SP+1) = 2, ..., P((N-1)*SP+1) = N
C     AND THEN REARRANGES P(1), P(SP+1), ..., P((N-1)*SP+1) SO THAT
C     A( (P(I)-1)*SA+1 ) .LE. A( (P(J)-1)*SA+1 ) IF AND ONLY IF
C     I .LT. J, WHERE I AND J SUBSCRIPT PROPER ELEMENTS OF P
C
      INTEGER SP, P(SP, *), SA, H, PIH, PI
      INTEGER A(SA, *)
C
C     CHECK INPUT PARAMETERS AND INITIALIZE H
C
      CALL I1SRT( SA, SP, N )
      IF ( I0SRT( 1, N, H ) .LT. 1 ) RETURN
C
C     INITIALIZE P
C
      DO 10 I = 1, N
 10     P(1, I) = I
C
C       CHECK IF DONE WITH SORT
C
 20     IF ( H .LT. 1 ) RETURN
        K = N - H
C
C       COMPARE
C
        DO 40 J = 1, K
          I = J
 30         IH = I + H
            PI = P(1, I)
            PIH = P(1, IH)
            IF ( A(1, PI) .LE. A(1, PIH) ) GOTO 40
C
C           EXCHANGE
C
            P(1, I) = PIH
            P(1, IH) = PI
C
C           PERCOLATE EXCHANGED LIST ELEMENT UP TO PROPER PLACE
C
            I = I - H
            IF ( I .GE. 1 ) GOTO 30
 40       CONTINUE
C
        H = ( H - 1 ) / 3
        GOTO 20
C
      END    
      SUBROUTINE I1SRT ( SA, SP, N )
C
C     I1SRT CHECKS LEGALITY OF VALUES OF SA, SP, N
C
      INTEGER SA, SP
C
C/6S
C     IF ( N .LT. 0 )
C    1   CALL SETERR( 27HSRTXXX - ILLEGAL VALUE OF N, 27, 1, 2 )
C     IF ( SA .LE. 0 )
C    1   CALL SETERR( 28HSRTXXX - ILLEGAL VALUE OF SA, 28, 2, 2 )
C     IF ( SP .LE. 0 )
C    1   CALL SETERR( 28HSRTXXX - ILLEGAL VALUE OF SP, 28, 3, 2 )
C/7S
C     IF ( N .LT. 0 )
C    1   CALL SETERR( 'SRTXXX - ILLEGAL VALUE OF N', 27, 1, 2 )
C     IF ( SA .LE. 0 )
C    1   CALL SETERR( 'SRTXXX - ILLEGAL VALUE OF SA', 28, 2, 2 )
C     IF ( SP .LE. 0 )
C    1   CALL SETERR( 'SRTXXX - ILLEGAL VALUE OF SP', 28, 3, 2 )
C/
C
      RETURN
C
      END


      INTEGER FUNCTION I0SRT ( SA, N, H )
C
C     I0SRT CHECKS INPUT PARAMETERS N, SA AND CALCULATES H
C     RETURNS H = 0 IF NO SORTING NECESSARY, ELSE
C     RETURNS SPACING, H, FOR FIRST INSERTION SORT.
C     I0SRT RETURNS TOTAL NUMBER OF ELEMENTS IN ARRAY = N * SA
C
      INTEGER SA, H
C
C/6S
C     IF ( N .LT. 0 )
C    1   CALL SETERR( 27HSRTXXX - ILLEGAL VALUE OF N, 27, 1, 2 )
C     IF ( SA .LE. 0 )
C    1   CALL SETERR( 28HSRTXXX - ILLEGAL VALUE OF SA, 28, 2, 2 )
C/7S  
C     IF ( N .LT. 0 )
C    1   CALL SETERR( 'SRTXXX - ILLEGAL VALUE OF N', 27, 1, 2 )
C     IF ( SA .LE. 0 )
C    1   CALL SETERR( 'SRTXXX - ILLEGAL VALUE OF SA', 28, 2, 2 )
C/
C
C     CHECK IF SORTING IS NECESSARY
C
      H = 0
      I0SRT = N * SA
      IF ( N .LE. 1 ) RETURN
C
C     CALCULATE H USING H NEW = 3 * H OLD + SA
C
      H = 4 * SA
C
 10     H = 3 * H + SA
        IF ( H .LT. I0SRT ) GOTO 10
C
      H = ( H - 4 * SA ) / 9
C
      RETURN
C
      END

