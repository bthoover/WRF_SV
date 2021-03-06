C**********************************************************************
C**********************************************************************
C
C
C
C	BELOW FOLLOW FOUR ROUTINES TO FACILITATE USE
C	OF THE LANCZOS PACKAGE NAMED LANCZ3.F
C
C	THIS IS NOW INTENDED TO SUPPLEMENT THE PACKAGE LANCZ4.F, 
C	AND THE CHANGE HAS BEEN MADE IN VRETR TO ALSO RETRIEVE
C	THE EIGENVALUES AND BOUNDS TO ALLOW DOUBLE CHECKING. 6-6-93.
C
C	MARTIN EHRENDORFER, NCAR, 18 MAY 1993
C
C	NAME OF ROUTINES:
C
C	CLANC
C	VRETR
C	OP
C	OPM
C	EVCHK   (ADDED ON 22 SEPTEMBER 1993) 
C
C
C
C**********************************************************************
	SUBROUTINE CLANC (N,LANMAX,LUNIT,KAPPA,W,NW,
     1             LSTEPS,NEIG,EV,BNDEV,NVEC,INDEV,INDVEC,IFAIL)
C**********************************************************************
C
C	SUBROUTINE TO CALL THE LANCZOS DRIVER LANDR. 
C	MARTIN EHRENDORFER, NCAR, 11 MAY 1993. 
C
C	THIS SUBROUTINE IS DEVISED TO MAKE MAXIMUM AND SIMPLE 
C	USE OF LANDR WITH A MINIMUM OF CHANGES IN THE 
C	LANDR-PACKAGE. FOR THE CHANGES IN THE LANDR-PACKAGE SEE
C	COMMENTS THERE. THESE CHANGES AFFECTED MAINLY THE 
C	COMPUTATION OF THE MACHINE EPSILON AND THE RANDOM
C	NUMBER GENERATOR.
C
C	INPUT:
C=====================================
C
C	N      ... SIZE OF THE PROBLEM
C	LANMAX ... MAXIMUM NUMBER OF LANCZOS STEPS ALLOWED
C	LUNIT  ... EIGENVECTOR INDICATOR, I/O UNIT
C	KAPPA  ... ACCURACY TO BE USED FOR VECTORS (IF LUNIT.GT.ZERO).
C	W      ... WORK ARRAY OF SIZE NW.
C       NW     ... AT LEAST 6*N+1+4*LANMAX IF LUNIT.LE.0
C	           AT LEAST 6*N+1+4*LANMAX+LANMAX*LANMAX IF LUNIT.GT.0
C
C	OUTPUT:
C=====================================
C	
C	KAPPA  ... ACCURACY ACTUALLY USED FOR EIGENVECTORS 
C	LSTEPS ... NUMBER OF LANCZOS STEPS ACTUALLY TAKEN
C	NEIG   ... NUMBER OF EIGENVALUES STABILIZED. THEY ARE TAKEN
C                  AS STABILIZED IF BNDEV(I).LE.16*EPS*EV(I). 
C	EV     ... EIGENVALUES OF THE PROBLEM. 
C		   THIS IS A VECTOR OF LENGTH LANMAX, CONTAINING 
C	           THE EIGENVALUES IN DECREASING ORDER. ONLY (THE
C		   FIRST) NEIG OF THEM ARE STABILIZED (SEE BELOW).
C	BNDEV  ... BOUNDS ON THE EIGENVALUES. VECTOR LIKE EV.
C	NVEC   ... NUMBER OF E-VECTORS WRITTEN ON LUNIT, NVEC.GE.NEIG.
C		   THEY ARE WRITTEN IF: BNDEV(K).LE.KAPPA*EV(K).
C	INDEV  ... INDICATOR FOR EIGENVALUES.
C		   THIS IS A VECTOR OF LENGTH LANMAX, SET TO ZEROES
C		   AND ONES AFTER SUCCESSFUL EXIT. A ONE IN POSITION
C	   	   I INDICATES ACCEPTANCE OF THE I-TH EIGENVALUE.
C		   IT WILL CONTAIN NEIG ONES. HOWEVER, IT IS 
C		   CONCEIVABLE THAT THE ONES DO NOT OCCUR IN THE 
C		   FIRST NEIG POSITIONS.
C	INDVEC ... INDICATOR FOR EIGENVECTORS.
C		   THIS IS A VECTOR OF LENGTH LANMAX, SET TO ZEROES
C		   AND ONES AFTER SUCCESSFUL EXIT IF LUNIT.GT.0. 
C		   IF LUNIT.LE.0 IT WILL CONTAIN ONLY ZEROES. 
C		   A ONE IN POSITION I SIGNIFIES THAT THE 
C		   VECTOR CORRESPONDING TO THE I-TH EIGENVALUE
C		   HAS BEEN WRITTEN TO LUNIT. OTHERWISE NO WRITING
C		   OCCURRED FOR THIS EIGENVALUE. 
C		   THE REASON FOR THESE INDEX CONTROLS IS: FIRST
C		   IT IS CONCEIVABLE THAT NOT THE FIRST NEIG 
C		   EIGENVALUES HAVE BEEN ACCEPTED, BUT SEVERAL MIGHT
C		   HAVE BEEN SKIPPED (THIS IS INDICATED BY A ZERO
C		   IN INDEV). SECOND, NOT THE FIRST NVEC VECTORS
C		   MIGHT HAVE BEEN WRITTEN, BUT SKIPPING MIGHT 
C		   HAVE OCCURRED. THIRD, DIFFERENCES MIGHT EXIST
C		   IN INDEXED POSITIONS. 
C		   SO, WHEN ASSOCIATING VECTORS AND VALUES, THESE
C 		   INDICES MUST BE CAREFULLY CHECKED. THIS IS NOT DONE
C		   HERE TO KEEP GENERALITY. 
C	IFAIL  ... ERROR FLAG. UNLESS IFAIL = 0, AN ERROR OCCURRED. 
C
	INTEGER NVEC,I
	INTEGER N,LANMAX,MAXPRS,LUNIT,LSTEPS,NEIG,NW,IFAIL,MSGLVL
	DOUBLE PRECISION CONDM,ENDL,ENDR,KAPPA,EPS
	DOUBLE PRECISION EV (LANMAX), BNDEV (LANMAX) , W (NW)
	INTEGER INDEV (LANMAX), INDVEC (LANMAX)
C
C	MORE SPECIFIC COMMENTS: 
C=====================================
C
C	ALWAYS USE LANMAX LESS OR EQUAL N.
C
C	IF LUNIT.LE.0 NO EIGENVECTORS ARE COMPUTED. OTHERWISE, LUNIT
C	IS TAKEN AS A NUMBER OF AN UNFORMATTED SEQUENTIAL ACCESS
C	FILE, TO WHICH THE EIGENVECTORS ASSOCIATED WITH CERTAIN 
C	EIGENVALUES ARE WRITTEN. THEY ARE WRITTEN IF: 
C	BNDEV(K).LE.KAPPA*EV(K). USE ROUTINE VRETR TO RETRIEVE
C	SELECTED EIGENVECTORS OF THOSE WRITTEN.
C
C	LANDR RETURNS THE EIGENVALUES IN ASCENDING ORDER. THIS 
C	ORDER IS REVERSED IN THIS ROUTINE AFTER SUCCESSFUL 
C	EXIT FROM LANDR. SO CLANC RETURNS THE EIGENVALUES OF
C	THE PROBLEM SORTED FROM LARGEST TO SMALLEST IN THE 
C	VECTOR EV. 
C
C	SPECIFICATION OF THE OPERATOR:
C=====================================
C
C	WHEN CALLING CLANC/VRETR TO SOLVE A STANDARD EIGENPROBLEM,
C	THE USER MUST SUPPLY A ROUTINE NAMED ATX93 ( X , N ) WHICH 
C	REPLACES THE N-DIMENSIONAL INPUT VECTOR X BY THE MATRIX
C	PRODUCT A TIMES X, WHERE A IS THE POSITIVE-DEFINITE, 
C	SYMMETRIC MATRIX/OPERATOR DEFINING THE EIGEN PROBLEM. 
C
C	NOTE THAT THE NAME OF THIS ROUTINE IS FIXED. 
C
C	IT IS IMPORTANT THAT THE OPERATOR A IS (POSITIVE 
C	SEMI-DEFINITE AND) SYMMETRIC. BOTH CONDITIONS ARE SATISFIED
C	IF A IS OBTAINED AS THE PRODUCT B(BT) OR (BT)B. THESE
C	CONDITIONS ARE NOT CHECKED BY THE LANCZOS ROUTINES, BUT 
C	ARE ASSUMED TO BE SATISFIED. 
C
C	THE ARGUMENTS OF ATX93 ( X , N ) ARE: 
C
C	X  ...  DOUBLE PRECISION VECTOR OF LENGTH N. ON INPUT 
C               THIS IS THE VECTOR ON WHICH THE MATRIX HAS TO OPERATE. 
C		ON OUTPUT, X MUST BE REPLACED BY 
C	        THE PRODUCT OF A TIMES X. 
C
C	N  ...  DIMENSION OF THE VECTOR X. MUST NOT BE CHANGED
C		WITHIN ATX. 
C
C	WITHIN ATX93 ONE MIGHT WANT TO WISH TO CHECK THE ADJOINT.
C	THIS IS EASY: (X , A(TRANS) A X) = (A X, A X)
C	NOTE : A(TRANS) A X IS THE RESULT OF THIS ROUTINE
C	AX IS THE RESULT OF THE OPERATOR A (E.G., TLM) ON X. 
C	A(TRANS) AX IS THE RESULT OF TRANSPOSED A (E.G., THE 
C	ADJOINT OF THE TLM) ON AX. THIS TEST IS NECESSARY AND 
C	(ALMOST) SUFFICIENT. 
C
C	EXAMPLE FOR ATX93. THIS EXAMPLE IMPLEMENTS ONE TIME STEP
C	OF THE FORWARD TIME CENTERED SPACE (FTCS) SCHEME FOR THE
C	LINEAR ADVECTIVE EQUATION. THIS SCHEME IS OF COURSE
C	UNCONDITIONALLY UNSTABLE WHICH SHOULD BE REFLECTED BY 
C	THE EIGENVALUES OF A. INDEED, THE EIGENVALUES OF A FOR
C	N=8 ARE: 4.532089, 4.532089, 3.347296, 3.347296, 2.000000, 
C	2.000000, 1.120615, 1.120615. 
C
C	SUBROUTINE ATX93 (X,N)
C	DOUBLE PRECISION X (*)
C	DOUBLE PRECISION Y (1000) .... THESE ARE WORKING ARRAYS
C	CALL  FWD (X,Y,N)         .... THIS PUTS AX ON Y
C	CALL AFWD (X,Y,N)         .... THIS PUTS ATRANS Y ON X
C	RETURN                    .... X IS REPLACED BY AT A X
C	END
C
C	SUBROUTINE FWD (X,Y,N)    .... Y IS THE SOLUTION AT THE 
C	DOUBLE PRECISION X(*), Y(*)    NEXT TIME STEP. ZERO VALUES
C	Y (1) = X (1) + X (2)          AT BOUNDARIES ASSUMED. 
C	Y (N) = X (N) - X (N-1)
C	DO 10 I = 2 , N-1
C10	Y (I) = X (I) + X (I+1) - X (I-1)
C	RETURN
C	END
C
C	SUBROUTINE AFWD (X,Y,N)   .... THIS IS THE ADJOINT OF FWD
C	DOUBLE PRECISION X(*), Y(*)
C	X (1) = Y (1) - Y (2) 
C	X (N) = Y (N) + Y (N-1)
C	DO 10 I = 2 , N-1
C10	X (I) = Y (I) - Y (I+1) + Y (I-1)
C	RETURN
C	END
C
C
C	SPECIAL SETTINGS IN CLANC: 
C=====================================
C
C	MAXPRS IS SET TO LANMAX, BECAUSE IT DOES NOT REALLY 
C	AFFECT THE COURSE OF ACTION IN LANDR. ONLY, IF MAXPRS
C	IS TOO SMALL, IT MAY HAPPEN THAT AN EXIT OCCURS BEFORE
C	THE NUMBER OF LANCZOS STEPS ONE IS WILLING TO DO HAS
C	BEEN DONE. 
C
	MAXPRS = LANMAX
C
C	ASSUMING A STANDARD EIGENVALUE SOLUTION IS SOUGHT, SET
C	CONDM (OF THE MATRIX M, THIS IS IDENTITY) TO ONE.
C
	CONDM = 1.0D00
C
C	SET THE INTERVAL TO SMALL NUMBERS CLOSE TO ZERO. THE 
C	ONLY WAY THAT THIS INTERVAL AFFECTS THE SOLUTION IS 
C	IN LANDR. NAMELY, AN EXIT FROM LANDR OCCURS WHENEVER
C	A STABILIZED EIGENVALUE IS WITHIN THE INTERVAL. THIS 
C	IS NOT DESIRABLE, IF ONE IS WILLING TO DO LANMAX 
C	ITERATIONS. THEREFORE, MAKE THIS INTERVAL SMALL. 
C	IN THIS CONTEXT THERE IS AN INCONSISTENCY IN THE 
C	DESCRIPTION OF LANDR, INSOFAR AS AN EXIT OCCURS
C	WHENEVER AN ACCEPTED EIGENVALUE IS WITHIN (ENDL,ENDR).
C
	ENDL = 0.0D00
	ENDR = 0.1D-04
C
C	KAPPA IS AN INPUT/OUTPUT ARGUMENT. IT IS ONLY USED/CHANGED
C	IF LUNIT.GT.0. ITS INPUT VALUE IS AFFECTED BY THE
C	STATEMENT IN LANDR: KAPPA = MAX(KAPPA,EPS**3/4). THUS, A
C	LARGER VALUE OF KAPPA WILL RESULT IN MORE EIGENVECTORS
C	BEING WRITTEN TO LUNIT IN RITVEC. IT IS RECOMMENDED TO TAKE A
C	RATHER SMALL VALUE OF KAPPA (E.G., 1.0D-30) TO BE 
C	CONSISTENT WITH THE CRITERION ON THE EIGENVALUES 
C	IN LANSO. 
C	NOTE THAT EPS = 0.25243548967D-28 ON THE CRAY-YMP. 
C	KAPPA IS USED IN RITVEC TO DETERMINE WHICH EIGENVECTORS 
C	ARE WRITTEN OUT. SINCE EPS**3/4 IS ALWAYS BIGGER THAN 
C	16*EPS, THE NUMBER OF EIGENVECTORS WRITTEN (NVEC) 
C	IS ALWAYS LARGER THAN THE NUMBER OF STABILIZED 
C	EIGENVALUES (NEIG).
C	MAKE KAPPA LARGER TO GET MORE EIGENVECTORS. 
C
C	MAKE CALL TO LANCZOS PACKAGE
C=====================================
C
	CALL LANDR ( N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,LUNIT,KAPPA,
     1               LSTEPS,NEIG,EV,BNDEV,W,NW,IFAIL,MSGLVL )
C
C	ON SUCCESSFUL EXIT FROM LANDR DO THREE THINGS: 
C=====================================
C
C	(1) CHANGE ORDER OF EIGENVALUES AND BOUNDS FROM LARGEST
C	TO SMALLEST. 
C
C	(2) DETERMINE THE NUMBER OF EIGENVECTORS WRITTEN OUT.
C	THIS COULD BE DONE MORE DIRECTLY 
C	IN RITVEC BUT IS DONE HERE FOR 
C	MINIMAL CHANGES IN THE LANDR-PACKAGE. THIS INFORMATION 
C	IS NEEDED FOR ROUTINE VRETR.
C
C	(3) TAKE CARE OF THE INDEX VECTORS INDEV AND INDVEC. 
C
	IF (IFAIL.EQ.0) THEN 
C
	DO 10 I = 1 , LANMAX
	W (I) = EV (I)
	W (I+LANMAX) = BNDEV (I)
10	CONTINUE
C
	DO 20 I = 1 , LANMAX
	EV (I) = W (LANMAX-I+1)
	BNDEV (I) = W (2*LANMAX-I+1)
20	CONTINUE
C
	DO 30 I = 1 , LANMAX 
	INDEV  (I) = 0
	INDVEC (I) = 0 
30	CONTINUE
C
	NVEC = 0 
	IF (LUNIT.GT.0) THEN 
	DO 40 I = 1 , LANMAX
	IF (BNDEV(I).LE.KAPPA*ABS(EV(I))) THEN 
	NVEC = NVEC + 1
	INDVEC (I) = 1
	ENDIF
40	CONTINUE
	ENDIF
C
	CALL DETEPS (EPS)
	DO 50 I = 1 , LANMAX
	IF (BNDEV(I).LE.16.0D00*EPS*ABS(EV(I))) INDEV (I) = 1
50	CONTINUE
C
	ENDIF
C
	RETURN
	END
C
C
C
C
C
C
C**********************************************************************
	SUBROUTINE VRETR ( L1,L2,LUNIT,NVEC,N,
     1                     EV,BNDEV,EVEC,IEVEC,JEVEC,IFAIL )
C**********************************************************************
C
C	SUBROUTINE TO RETRIEVE (SOME OF) THE EIGENVECTORS
C	AFTER A SUCCESSFUL CALL TO CLANC. 
C	MARTIN EHRENDORFER, NCAR, 11 MAY 1993. 
C
C	AFTER A SUCCESSFUL CALL TO CLANC, THE UNIT LUNIT
C	CONTAINS NVEC VECTORS. THE FIRST OF THESE IS THE EIGENVECTOR
C	ASSOCIATED WITH THE SMALLEST ACCEPTED EIGENVALUE (ACCEPTED
C	FOR COMPUTATION OF THE VECTOR), AND THE LAST (I.E., AT THE NVEC 
C	POSITION) IS THE EIGENVECTOR ASSOCIATED WITH THE LARGEST 
C	EIGENVALUE. THIS WRITING IS DONE IN RITVEC. 
C
C	SINCE WE ARE INTERESTED IN THE LARGEST EIGENVECTORS, THIS 
C	ROUTINE DOES IMPLICITLY SOME REORDERING. 
C
C	THIS ROUTINE RETRIEVES FROM LUNIT ALL
C	EIGENVALUES WITH INDEX L = L1 , L2, 1
C	WHERE L1 .LE. L2, AND THE EIGENVALUE OF INDEX L1 IS 
C	LARGER THAN THE EIGENVALUE ASSOCIATED WITH INDEX L2. 
C	FOR EXAMPLE, L1=2, L2=5, SIGNIFIES TO THIS ROUTINE THAT
C	WE WANT TO RETRIEVE THE FOUR EIGENVECTORS ASSOCIATED 
C	WITH THE SECOND LARGEST EIGENVALUE UP TO THE FIFTH 
C	LARGEST EIGENVALUE. 
C
C	ON SUCCESSFUL EXIT (IFAIL = 0), 
C	THE EIGENVECTORS ARE RETURNED AS COLUMNS OF EVEC. THE 
C	FIRST COLUMN OF EVEC CONTAINS THE EIGENVECTOR ASSOCIATED
C	WITH THE EIGENVALUE OF INDEX L1, AND THE LAST COLUMN
C	THE VECTOR ASSOCIATED WITH INDEX L2. THUS, THE COLUMNS
C	CONTAIN THE VECTORS CORRESPONDING TO DECREASING EIGENVALUES
C	FROM L1 TO L2. SPECIFICALLY, IF L1=1 AND L2=NVEC, THEN 
C	EVEC WILL CONTAIN ALL ACCEPTED EIGENVECTORS, FROM LARGEST
C	TO SMALLEST, IN COLUMNS 1 TO L2-L1+1=NVEC. 
C
C	CHANGE ON 6-6-93: WHEN READING THE VECTORS FROM LUNIT
C	NOW ALSO THE EIGENVALUES PLUS BOUNDS ARE READ AND RETURNED
C	TO THE CALLING ROUTINE. THIS ALLOWS FOR DOUBLE CHECKING
C	OF THE VECTOR/VALUE RELATIONSHIP IN THE CALLING ROUTINE. 
C	DONE 6-6-93. 
C
	INTEGER L1, L2, LUNIT, NVEC, N, IEVEC, JEVEC, IFAIL
	DOUBLE PRECISION EVEC (IEVEC,JEVEC) , EV (*), BNDEV (*)
C
C	INPUT: 
C=====================================
C
C	L1    ... INDEX OF FIRST EIGENVECTOR DESIRED
C	L2    ... INDEX OF LAST  EIGENVECTOR DESIRED 
C		  L2.GE.L1 AND L2.LE.NVEC
C	LUNIT ... UNIT NUMBER USED FOR CLANC (LUNIT.GT.0)
C	NVEC  ... NUMBER OF EIGENVECTORS WRITTEN AS RETURNED BY CLANC
C	N     ... DIMENSION OF THE PROBLEM
C	IEVEC ... FIRST DIMENSION OF EVEC AS IN CALLER
C	JEVEC ... SECOND DIMENSION OF EVEC AS IN CALLER
C
C	OUTPUT:
C=====================================
C
C	EVEC  ... COLUMNS 1 TO L2-L1+1 CONTAIN THE DESIRED E-VECTORS.
C	EV    ... ARRAY CONTAINING E-VALUES FOR COLUMNS OF EVEC.
C	BNDEV ... BOUNDS FOR E-VALUES IN EV.
C	IFAIL ... ERROR FLAG, UNLESS IFAIL=0, AN ERROR OCCURRED. 
C
C
	INTEGER IDUMY,JDUMY,M2,I,J,M1
	DOUBLE PRECISION XDUMY,YDUMY
C
C	CHECK FOR INPUT CONSISTENCY
C=====================================
C
	M2 = L2 - L1 + 1
	IFAIL = -1 
	IF (L1.GT.L2)           RETURN
	IF (LUNIT.LE.0)         RETURN
	IF (NVEC.EQ.0)          RETURN
	IF (IEVEC.LT.N)         RETURN
	IF (JEVEC.LT.M2)        RETURN
	IFAIL = 0 
C
C	READ THROUGH THE FILE AND STORE THE REQUIRED EIGENVECTORS
C=====================================
C
	OPEN (LUNIT,FORM='UNFORMATTED',STATUS= 'UNKNOWN')
	REWIND (LUNIT)
	READ (LUNIT) IDUMY, JDUMY, XDUMY
	DO 10 I = 1 , NVEC
	J = NVEC - I + 1
	READ (LUNIT) EV(M2) , BNDEV(M2) , (EVEC(M1,M2),M1=1,N)
	IF (J.LE.L2) M2 = M2 - 1
	IF (M2.EQ.0) GOTO 20
10	CONTINUE
20	CLOSE (LUNIT)
C
	RETURN
	END
C
C
C
C
C
C
C**********************************************************************
	SUBROUTINE OP (N,X,Y,Z)
C**********************************************************************
C
C	SUBROUTINE TO DEFINE THE OPERATOR FOR WHICH THE
C	SPECTRAL DECOMPOSITION IS DESIRED. 
C	MARTIN EHRENDORFER, NCAR, 14 MAY 1993.
C
C	THIS IS THE ROUTINE ORIGINALLY REQUIRED BY THE 
C	LANCZOS PACKAGE. IN THIS IMPLEMENTATION THE USER
C	DOES NO LONGER SUPPLY OP, BUT ATX93 TO DEFINE THE 
C	OPERATOR A. 
C
C	IN ORDER TO SOLVE THE GENERALIZED EIGENPROBLEM: 
C
C	K X = LAMBDA M X
C
C	THIS ROUTINE HAS TO SOLVE THE EQUATION:
C
C	M R = K Q 
C
C	FOR GIVEN Q. FOR A STANDARD PROBLEM M IS IDENTITY AND THE
C	RESULT OF THIS ROUTINE IS SIMPLY THE OPERATOR K APPLIED ON Q. 
C
C	THIS STANDARD PROBLEM IS SOLVED IN THE SITUATION CONSIDERED
C	HERE / BY THE PROGRAMS DESCRIBED HERE. 
C
C	IN THIS ROUTINE Z IS THE OUTPUT, AND X AND Y ARE
C	THE INPUT. IN ALL CASES I RAN, X AND Y ARE THE SAME,
C	SO IT DOES NOT REALLY MATTER, ONTO WHICH TO APPLY 
C	THE OPERATOR. FOR THE STANDARD PROBLEM I SUGGEST TO SET:
C
C	Z = K X
C
C	WHERE X IS INPUT VECTOR, K IS THE OPERATOR, Z IS RESULT.
C
	DOUBLE PRECISION X (*) , Y (*) , Z (*)
C
C	PUT THE INPUT VECTOR X ONTO Z
C
	DO 10 I = 1 , N
10	Z (I) = X (I) 
C
C	HAVING DONE THIS, REPLACE Z BY A Z ( THAT IS K Z )
C	THROUGH A CALL TO THE USER SUPPLIED ROUTINE ATX93
C
	CALL ATX93 ( Z , N )
C
C	RETURN TO THE LANCZOS ROUTINES 
C
	RETURN
	END
C
C
C
C
C
C**********************************************************************
	SUBROUTINE OPM (N,X,Y)
C**********************************************************************
C
C	SUBROUTINE TO RETURN Y = M X FOR GIVEN X. 
C	IN THE STANDARD PROBLEM M IS IDENTITY. 
C	MARTIN EHRENDORFER, NCAR, 18 MAY 1993.
C
C	IN THIS ROUTINE X IS THE INPUT AND Y IS THE OUTPUT. 
C	Y IS TO BE SET TO M X.
C
C	SINCE WE ARE DEALING WITH THE STANDARD EIGENPROBLEM
C	Y IS SIMPLY REPLACED BY X. THE USER SHOULD NOT MODIFY 
C	THIS ROUTINE. 
C
	DOUBLE PRECISION X (*) , Y (*)
C
	DO 10 I = 1 , N 
10	Y (I) = X (I) 
C
	RETURN
	END
C
C
C
C
C**********************************************************************
C**********************************************************************
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	SUBROUTINE EVCHK ( N,M,LUNIT1,LUNIT2,ACC,
     1                     X1,X2,X3,X4,W,NW,IFAIL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	THIS ROUTINE IS DESIGNED TO CHECK THE EIGENVALUE/VECTOR
C	PAIRS THAT ARE PRODUCED BY THE LANCZOS ALGORITHM.
C
C	MARTIN EHRENDORFER, NCAR, 22 SEPTEMBER 1993. 
C
C	CALLING THIS ROUTINE REQUIRES THE ROUTINE ATX93 THAT HAS BEEN 
C	USED TO PRODUCE THE EIGENPAIRS AS WELL AS THE DIRECT ACCESS
C	FILE ON WHICH THE EIGENPAIRS ARE STORED. 
C	ALSO THIS ROUTINE NEEDS THE ROUTINE DDOT. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DESCRIPTION OF TESTING:	
C
C	(1): CHECK ORTHOGONALITY OF EIGENVECTORS.
C
C	(2): SEND EIGENVECTORS THROUGH THE MATRIX.
C
C	BOTH TESTS ARE PERFORMED IF M.GT.1. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	THIS ROUTINE SHOULD BE CALLED IN A NEW RUN AFTER THE 
C	LANCZOS ITERATION HAS BEEN SUCCESSFULLY COMPLETED. 
C	IN THIS CASE, LUNIT1 AND LUNIT2 MUST BE OPENED IN THE CALLING 
C	ROUTINE, AND LUNIT1 IS RESERVED FOR OUTPUT, AND LUNIT2 CONTAINS
C	THE EIGENVALUE/VECTOR PAIRS. SEE BELOW. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DESCRIPTION OF ARGUMENTS OF THIS ROUTINE:
C
C	INPUT: 
C
C	N           ...   SIZE OF THE PROBLEM (I.E. LENGTH OF VECTORS) 
C
C	M           ...   NUMBER OF E-VECTORS TO BE TESTED 
C
C       LUNIT1      ...   UNIT ONTO WHICH THE WRITTEN OUTPUT
C		          OF THIS ROUTINE IS MADE. IF LUNIT1 IS LESS OR
C			  EQUAL TO ZERO THEN NO WRITTEN OUTPUT IS MADE
C	  		  IN THIS ROUTINE. LUNIT1 MUST BE OPENED
C			  IN CALLING ROUTINE. 
C
C	LUNIT2      ...   UNIT OF DIRECT ACCESS FILE CONTAINING 
C		          THE EIGENVECTORS. LUNIT2 MUST BE OPENED
C			  IN CALLING ROUTINE. IT IS ASSUMED THAT THE 
C			  J-TH EIGENVALUE/VECTOR PAIR HAS BEEN WRITTEN 
C			  TO THIS FILE WITH THE STATEMENT: 
C			  WRITE (LUNIT2,REC=J) EV(J), (X(I),I=1,N) 
C
C	ACC 	    ...   ACCURACY CRITERION, SET TO A SMALL NUMBER. 
C			  
C	X1,X2,X3,X4 ...   WORKING VECTORS OF LENGTH N
C
C	W           ...   WORKING ARRAY OF LENGTH NW
C
C	NW	    ...   LENGTH OF W, AT LEAST 6*N
C
C
C	OUTPUT: 
C
C	IFAIL       ...   ERROR INDICATOR, WILL BE ZERO IF 
C			  ALL THE DESIRED TESTS ARE SATISFIED
C
C	NOTE: IF LUNIT IS GT ZERO WRITTEN OUTPUT IS MADE TO LUNIT1
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	REAL X1 (N), X2 (N), X3 (N), X4 (N), W (NW)
C
	IFAIL = 0
	IF ( M.LT.1 ) RETURN
	IFAIL = 0
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	LEVEL 1 TESTING (ORTHOGONALITY OF EIGENVECTORS) 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	IF ( LUNIT1 .GT. 0 ) THEN 
        WRITE (LUNIT1,'(/A/)') 'RESULTS FROM ORTHOGONALITY TESTING'
	ENDIF
C
	DO 20 J = 1 , M
C
C	READ ONE EIGENVECTOR / EIGENVALUE PAIR ONTO X1 (I.E. XJ)
C
	READ (LUNIT2,REC=J) W (1) , W (2) , ( X1 (I) , I = 1 , N )
C
C	COMPUTE XK TRANS XJ AND CHECK THE RESULT
C
	DO 21 K = J , M
	READ (LUNIT2,REC=K) W (1) , W (2) , ( X2 (I) , I = 1 , N )
	X3 (K) = DDOT ( N,X2,1,X1,1 ) 
	A = 0.0
	IF ( K.EQ.J ) A = 1.0
	IF ( ABS (X3 (K)-A) .GT. ACC ) IFAIL = IFAIL + 1
21	CONTINUE
C
	IF ( LUNIT1 .GT. 0 ) THEN 
	WRITE ( LUNIT1,300 ) J , ( X3 (L) , L = J , M ) 
300	FORMAT (/2X,I7,5E15.7/40(9X,5E15.7/)/) 
	ENDIF
C
20	CONTINUE
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	LEVEL 2 TESTING (ORTHOGONALITY WITH RESPECT TO A) 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	IF ( LUNIT1.GT.0 ) THEN 
	WRITE (LUNIT1,'(/A/)') 'RESULTS FROM ORTHOGONALITY WRT A TEST'
	ENDIF
C
	DO 30 J = 1 , M
C
C	READ ONE EIGENVECTOR / EIGENVALUE PAIR ONTO X1 (I.E. XJ)
C
	READ (LUNIT2,REC=J) W (1) , W (2) , ( X1 (I) , I = 1 , N )
C
C	SAVE THE CURRENT EIGENVECTOR ON X2 AND VALUE ON XLMJ
C
	DO 31 I = 1 , N
31	X2 (I) = X1 (I) 
C
	XLMJ = W (1) 
C
C	SEND THE CURRENT EIGENVECTOR THROUGH THE MATRIX
C
	CALL ATX93 ( X1 , N ) 
C
C	NOW: X1 IS: A XJ, AND X2 IS: XJ, AND XLMJ IS: LAMBDA J
C
C	(1): COMPUTE A XJ - LAMBDA J XJ
C	SAVE THE DOT PRODUCT OF THE ABOVE QUANTITY ON W (J+2)
C
	DO 32 I = 1 , N
	X3 (I)  = X1 (I) - XLMJ * X2 (I)
32	X4 (I)  = X3 (I)
	W  (J+2)  = DDOT ( N,X3,1,X4,1 )
 	IF ( ABS (W(J+2)) .GT. ACC ) IFAIL = IFAIL + 1
C
C	(2): COMPUTE QUOTIENTS
C	XJ T A XJ / XJ T XJ - LAMBDA J
C	SAVE THE RESULT ON W AGAIN 
C
	W  (J+N+2)  = DDOT ( N,X2,1,X1,1 )
	DO 33 I = 1 , N
33	X3 (I)  = X2 (I)
	W  (J+N+2)  = W (J+N+2) / DDOT ( N,X3,1,X2,1 ) - XLMJ
	IF ( ABS (W(J+N+2)) .GT. ACC ) IFAIL = IFAIL + 1
C
C	(3): COMPUTE XK TRANS A XJ
C
	DO 34 K = 1 , M 
	READ (LUNIT2,REC=K) W (1) , W (2) , ( X3 (I) , I = 1 , N )
	X4 (K) = DDOT ( N,X3,1,X1,1 )
	A = 0.0
	IF ( K.EQ.J ) A = XLMJ 
	IF ( ABS ( X4(K)-A ) .GT. ACC ) IFAIL = IFAIL + 1
34 	CONTINUE
C
	IF ( LUNIT1 .GT. 0 ) THEN 
        WRITE (LUNIT1,'(/A/)') 'RESULTS FROM XK TRANS A XJ TEST'
	WRITE ( LUNIT1,300 ) J , ( X4 (L) , L = 1 , M ) 
	ENDIF
C
30	CONTINUE
C
	IF ( LUNIT1 .GT. 0 ) THEN 
        WRITE (LUNIT1,'(/A/)') 'RESULTS FROM A XJ - LAM J XJ TEST'
	WRITE ( LUNIT1,300 ) 1 , ( W (J+2)   , J = 1 , M ) 
        WRITE (LUNIT1,'(/A/)') 'RESULTS FROM QUOTIENT TEST'
	WRITE ( LUNIT1,300 ) 2 , ( W (J+N+2) , J = 1 , M ) 
	ENDIF
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	END OF TESTING
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
	RETURN
	END
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C

