cc34567
      program read_output
      parameter(N=3600520,NP2=N+2,NJ=6)
      double precision EV(N), BNDEV(N), V_V(N,NJ)

      MMATRIX = N
C
        OPEN (UNIT=23,FILE='vector.dat',FORM='UNFORMATTED',
     1                                        STATUS='old',           
     2        ACCESS='DIRECT',RECL=NP2*8)
        DO 30 J = 1 , NJ
        read (23,REC=J) EV(J),BNDEV(J),(V_V(I,J),I=1,MMATRIX)
	write(6,*) ev(j)
	write(30+j)v_v(:,j)
30	continue
	end

