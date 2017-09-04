stepsize_n(-1)=0.5_dp
    itermax = getinput%int("maximum_iteration_nbr", defaultvalue=huge(1), assert=">0")
    i=0
    f=0._dp
    fold=huge(1._dp)
    fmin=huge(1._dp)
    oldDeltaF=huge(1._dp)
    j=1
    call energy_and_gradient(f,df)
    df_prev = df
    open(12,file="output/iterate.dat")
    
    do while(i<itermax)
      print*,""
      print*
      print*
      print*,"ITERATION", j
      fold=f
      oldDeltaF=deltaF
 
      ich_continue=.true.
      k=0
      find_new_value = .false.
      stepsize=stepsize_n(i-1)*1.5_dp
      
      do while(ich_continue)
        solvent(1)%xi=solvent(1)%xi-stepsize*df_prev(:,:,:,:,1)
        call energy_and_gradient(f, df)
        solvent(1)%xi=solvent(1)%xi+stepsize*df_prev(:,:,:,:,1)

        if(f<fmin) then
          stepsize_giving_minimum_f=stepsize
          df_for_steepsize_giving_minimum_f = df
          fmin=f
          deltaF = abs(fmin-fold)/max(fmin,abs(fold),1._dp)
          find_new_value = .true.
          if( (deltaF < factr) ) then
            ich_continue=.true.
          else
            ich_continue=.false.
          end if
        else
          if( k .lt. n_try_max_by_iteration) then
            ich_continue=.true.
            stepsize=stepsize*0.5_dp
          else
            ich_continue=.false.
          end if
        end if
        stepsize_n(i)=stepsize_giving_minimum_f
        PRINT*,"AT ITERATION ",j,"BEST STEPSIZE TO DATE=",stepsize_giving_minimum_f,"F=",f,"FMIN=",fmin,"ACTUAL STEPSIZE", stepsize
        k=k+1
      end do
      PRINT*,"AT ITERATION ",j,"I WILL USE STEPSIZE",stepsize_giving_minimum_f
      PRINT*
      PRINT*
      if( .not. find_new_value ) then
          print*, "Minimization finished without converging! criteria=", oldDeltaF, " for ", factr
          exit
      end if      
      
      solvent(1)%xi=solvent(1)%xi-stepsize_giving_minimum_f*df_prev(:,:,:,:,1)
      df_prev = df_for_steepsize_giving_minimum_f
      i=i+1
      j=j+1
      write(12,*) i,fmin
      if( deltaF < factr ) exit
    end do
    close(12)
