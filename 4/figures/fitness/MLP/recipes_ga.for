c
c
c
c
c-------------------------------------------------------------
c>>>>>>>>>>>>>>>>>>>>>>>  Subroutines from Numerical Recipes:
c-------------------------------------------------------------

      FUNCTION RAN3(IDUM)
        IMPLICIT REAL*8(M)
        PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
c      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
	IFF=1
	MJ=MSEED-IABS(IDUM)
	MJ=MOD(MJ,MBIG)
	MA(55)=MJ
	MK=1
	DO 11 I=1,54
	  II=MOD(21*I,55)
	  MA(II)=MK
	  MK=MJ-MK
	  IF(MK.LT.MZ)MK=MK+MBIG
	  MJ=MA(II)
11      CONTINUE
	DO 13 K=1,4
	  DO 12 I=1,55
	    MA(I)=MA(I)-MA(1+MOD(I+30,55))
	    IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
	INEXT=0
	INEXTP=31
	IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
      
      
      SUBROUTINE SORT3(N,RA,RB,RC,WKSP,IWKSP)
      implicit double precision (a-h,o-z)
      DIMENSION RA(N),RB(N),RC(N),WKSP(N),IWKSP(N)
      CALL INDEXX(N,RA,IWKSP)
      DO 11 J=1,N
	WKSP(J)=RA(J)
11    CONTINUE
      DO 12 J=1,N
	RA(J)=WKSP(IWKSP(J))
12    CONTINUE
      DO 13 J=1,N
	WKSP(J)=RB(J)
13    CONTINUE
      DO 14 J=1,N
	RB(J)=WKSP(IWKSP(J))
14    CONTINUE
      DO 15 J=1,N
	WKSP(J)=RC(J)
15    CONTINUE
      DO 16 J=1,N
	RC(J)=WKSP(IWKSP(J))
16    CONTINUE
      RETURN
      END
      
      
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit double precision (a-h,o-z)
      DIMENSION ARRIN(1),INDX(1)
c
      if (n.eq.1) then
        indx(1)=1
        return
      endif
c
      DO 11 J=1,N
	INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
	IF(L.GT.1)THEN
	  L=L-1
	  INDXT=INDX(L)
	  Q=ARRIN(INDXT)
	ELSE
	  INDXT=INDX(IR)
	  Q=ARRIN(INDXT)
	  INDX(IR)=INDX(1)
	  IR=IR-1
	  IF(IR.EQ.1)THEN
	    INDX(1)=INDXT
	    RETURN
	  ENDIF
	ENDIF
	I=L
	J=L+L
20      IF(J.LE.IR)THEN
	  IF(J.LT.IR)THEN
	    IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
	  ENDIF
	  IF(Q.LT.ARRIN(INDX(J)))THEN
	    INDX(I)=INDX(J)
	    I=J
	    J=J+J
	  ELSE
	    J=IR+1
	  ENDIF
	GO TO 20
	ENDIF
	INDX(I)=INDXT
      GO TO 10
      END

