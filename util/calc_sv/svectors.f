	PROGRAM SVECTORS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             MAIN PROGRAM FOR THE CALCULATION OF SINGULAR             C        
C             VECTORS ON MM5.                                          C
C             THE PROCEDURE IS BASED ON THE USE OF THE LANCZOS METHOD  C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C*****************************************************************
        PARAMETER (MMATRIX = 3600520)
	PARAMETER ( LANMAX = 25)
        PARAMETER ( LUNIT1 = 21, LUNIT2 = 22, LUNIT3 = 23 ) 
	PARAMETER ( NW = 6*MMATRIX+1+4*LANMAX+LANMAX*LANMAX )
	PARAMETER ( IV = MMATRIX , JV = 75 ) 
	PARAMETER ( NP2 = MMATRIX + 2 ) 
C
	REAL W_W ( NW ) 
	REAL V_V ( IV , JV ) 
C234567
	REAL X1 (MMATRIX), X2 (MMATRIX), X3 (MMATRIX),
     $                   X4 (MMATRIX)
C
	REAL EV (MMATRIX), BNDEV (MMATRIX)
	INTEGER INDEV (LANMAX), INDVEC (LANMAX) 
C
        REAL IC(MMATRIX), XKAPPA
        LOGICAL RUNLIN
C
        print*, 'in prog. svectors/ dim. of the pb =',MMATRIX
C
        RUNLIN=.TRUE.! IF TRUE, LINEAR MODEL WILL BE INTEGRATED FORWARD
C                      WITH INDIVIDUAL SV'S AS INITIAL CONDITIONS.
cmcm        print*, 'CALL PRE_ATX93'
cmcm        CALL PRE_ATX93 
C
	OPEN ( LUNIT1 , FILE = 'lanc.dat',STATUS='NEW',FORM='FORMATTED') 
C	
	XKAPPA = 1.0E-30
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmcm        print*, 'In main program NDIMS=',NDIMS
cmcm        print*, 'In main program MIX,MJX,MKX,MIXNH,MJXNH,KXP1NH,MKXNH'
cmcm        print*,  MIX,MJX,MKX,MIXNH,MJXNH,KXP1NH,MKXNH
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C====================================================================
	CALL CLANC ( MMATRIX,LANMAX,LUNIT2,XKAPPA,W_W,NW,
     1       LSTEPS,NEIG,EV,BNDEV,NVEC,INDEV,INDVEC,IFAIL ) 
C====================================================================
C
      WRITE (LUNIT1,'(A,I10)') 'DIMENSION OF PROBLEM MMATRIX=', MMATRIX
      WRITE (LUNIT1,'(A,I10)')    'MAX. # ITS. ALL. LANMAX = ', LANMAX
      WRITE (LUNIT1,'(A,I10)')    'UNIT F. LANCZOS LUNIT2  = ', LUNIT2
      WRITE (LUNIT1,'(A,E22.15)') 'ACC. WRITE VECTOR KAPPA = ', XKAPPA
      WRITE (LUNIT1,'(A,I10)')    'STEPS ACT. TAKEN LSTEPS = ', LSTEPS
      WRITE (LUNIT1,'(A,I10)')    'E-VALUES ACCEPTED NEIG  = ', NEIG
      WRITE (LUNIT1,'(A,I10)')    'E-VECTORS ACCEPTED NVEC = ', NVEC
      WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C
	WRITE (LUNIT1,201) ( EV     (I) , I = 1 , LANMAX )
	WRITE (LUNIT1,202) ( BNDEV  (I) , I = 1 , LANMAX )
201     FORMAT (/2X,'E-VALUS',5E15.7/40(9X,5E15.7/))
202     FORMAT (/2X,'E-BNDS ',5E15.7/40(9X,5E15.7/))
C
	WRITE (LUNIT1,204) ( INDEV   (I) , I = 1 , LANMAX )
	WRITE (LUNIT1,205) ( INDVEC  (I) , I = 1 , LANMAX )
204     FORMAT (2X,'INDEV  ',5I15/40(9X,5I15/))
205     FORMAT (2X,'INDVEC ',5I15/40(9X,5I15/))
C
C
C====================================================================
C	IF CLANC SUCCESSFUL CALL VRETR TO RETRIEVE EIGEN V E C T O R S 
C====================================================================
C
C
	IF ( IFAIL.EQ.0 ) THEN
C	
	L1 = LANMAX-LSTEPS+1 
        L2=  NVEC
C:M.PONDECA/13MARCH1997: FOR LARGE LANMAX AN UNDESIRABLE EXIT CAN OCCUR
C BEFORE THE LANMAX NUMBER OF LANCZOS STEPS ONE IS WILLING TO PERFORM.
C THIS MAY NEED TO BE FIXED! WHEN SUCH AN EXIT OCCURS THAN L1 HAS
C TO BE AT LEAST (LANMAX-LSTEPS+1), OTHERWISE VRETR WILL CRASH!
C NOTE THAT SUCH AN EXIT IS ALSO REFLECTED IN THE FIELD EV(LANMAX).
C THE FIRST (LANMAX-LSTEPS) VALUES IN THIS FIELD WILL BE ZEROS! THEN
C COMES THE LARGEST EIGENVALUE OF THE PROBLEM.
C
C====================================================================
      CALL VRETR(L1,L2,LUNIT2,NVEC,MMATRIX,EV,BNDEV,V_V,IV,JV,IFAIL)
C====================================================================
C
	WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C
	DO 20 I = 1 , MMATRIX
20	WRITE (LUNIT1,'(2X,3F15.6)') (V_V(I,J),J=1,3)
C
C	WRITE OUT THE EIGENVECTORS TO DIRECT ACCESS 
C====================================================================
C
	OPEN (UNIT=LUNIT3,FILE='vector.dat',FORM='UNFORMATTED',
     1                                        STATUS='NEW',           
     2        ACCESS='DIRECT',RECL=NP2*8)
	DO 30 J = 1 , L2-L1+1
30      WRITE (LUNIT3,REC=J) EV(J),BNDEV(J),(V_V(I,J),I=1,MMATRIX)
C
C	TEST THE EIGENVECTORS THAT HAVE BEEN RETRIEVED
C====================================================================
C
	MTEST = 3!  NVEC
        IF (NVEC.LT.3) MTEST=NVEC
        ACC = 1.0E-10
C
C====================================================================
       CALL EVCHK (MMATRIX,MTEST,LUNIT1,LUNIT3,
     -              ACC,X1,X2,X3,X4,W_W,NW,IFAIL)
C====================================================================
C
	WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C	
	ENDIF
C*****************************************************************
C        RETRIEVE TRUE SINGULAR VECTORS OF THE SYSTEM
C        THROUGH A LINEAR TRANSFORMATION.
C        THESE WILL BE STORED IN FILE 'results.out' IN 
C        SUBROUTINE 'true_svs'.
C*****************************************************************
C
C       OPEN (13, FILE='results.out',FORM='FORMATTED',STATUS='NEW')
C       OPEN (14, FILE='results2.out',FORM='FORMATTED',STATUS='NEW')
C         OPEN (8, FILE='results.log',FORM='FORMATTED',STATUS='NEW')
C         print*, 'WRITING OUT THE EIGENVALUES'
C         print*, 'J,       BNDEV (J),  AMPL. FACT. ENERGY ERROR'
C         WRITE (8,*) 'J,       BNDEV (J),  AMPL. FACT. ENERGY ERROR'
C         print*, '----------------------------'
C         WRITE (8,*) '----------------------------'
C         WRITE (13,*) (EV(J),J=1,L2-L1+1)
C         WRITE (13,*) (BNDEV(J),J=1,L2-L1+1)
C         WRITE (14,*) (EV(J),J=1,L2-L1+1)
C         WRITE (14,*) (BNDEV(J),J=1,L2-L1+1)
C
C         DO J=1,  L2-L1+1
C            CALL TRUE_SVS(V_V(1,J),MMATRIX,J,ENERGY,RUNLIN) 
C            IF (RUNLIN) ERROR=ABS((EV(J)-ENERGY)/EV(J))
C           IF (.NOT.RUNLIN) ERROR=9999.999
C            WRITE (6,1640) J,BNDEV(J),EV(J),ENERGY,ERROR
C            WRITE (8,1640) J,BNDEV(J),EV(J),ENERGY,ERROR
C         ENDDO
C1640     FORMAT(I4,1X,E12.6,3(1X,E12.6))
C
C======================================================================
        print*, ' ALL CALCULATIONS  DONE'
C
9999    STOP
	END
C 
c        include 'grad.f_improved'
C        include 'grad_SVS.f_JACOBIAN_correct'
CCCCC        include 'grad_SVS.f'
