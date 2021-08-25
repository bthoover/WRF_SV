
      SUBROUTINE ATX93(XIC,N)

      DOUBLE PRECISION XIC(N)
      print*,'+++++CHECK PARAMETERS IN ATX93+++++'
      print*,'In ATX93 MMAT=',N
      if(N.ne.3600520) then
         print*, 'hmm ..... interesting'
         stop
      end if
      print*,'in atx93_mcm: >>>> ', XIC(1), XIC(2), XIC(N-1), XIC(N)
      print*,'+++++END CHECK PARAMETERS IN ATX93+++++'
      call execute_command_line("rm fort.16 fort.17", exitstat=i)
      print *, "Exit status of remove state vector was ", i

      if(i.eq.0) print*, 'files fort.16 and fort.17 were removed'

      open(16,file='fort.16',status='new',form='unformatted')
      write(16) XIC
      
      print *, "in progress"

      call system("python3 govern_wrf_tlm_adj.py")
      open(17,file='fort.17',status='old',form='unformatted')
      read(17) XIC

      print *, "PYTHON: DONE"
      

      RETURN
      END 



    
