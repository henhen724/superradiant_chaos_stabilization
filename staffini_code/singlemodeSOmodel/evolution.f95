MODULE evolution
  USE global, ONLY: varsize, decode_state, write_output, nmax, mmax, &
                  & convergence_threshold, CurrentState, Lam, Phi, &
                  & tinit, tend, tstep
  USE eom_so, ONLY: diffeqn
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  INTEGER, PARAMETER :: BackCheck = 40
  REAL (KIND=DP) :: StateStore(BackCheck, 1:2)
  
  ! State whether to dump info relating to convergence or not
  LOGICAL :: dbg_convergence
  
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: execute_single
  ! VARIABLES:
  ! SYNOPSIS:
  !   This routine simply performs the evolution for the first value
  !   of P passed to the program.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  SUBROUTINE execute_single()
    USE global, ONLY: Pump, OmegaTil
    INTEGER (KIND=DP) :: buffer=3000
    
    write(*,*) Pump, OmegaTil
    call next_steady(CurrentState,tinit,tend,tstep, buffer)

  END SUBROUTINE execute_single

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: next_steady
  ! VARIABLES:
  !           state  - encoded state (INOUT)
  !           xinit  - initial state to pass to NAG (INOUT)
  !           xend   - final state from NAG (INOUT)
  !           eom    - subroutine to use as eqn of motion (INTERFACE)
  !           step   - timestep
  !           buffer - time before code starts checking for 
  !                                               convergence (IN)
  ! SYNOPSIS:
  !   This evolves the state until a steady state is reached 
  !   (|ds_dt|< tolerance), or returns an error if no steady
  !   state is reached by xend.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  SUBROUTINE next_steady(state, xinit, xend, step, buffer)

    REAL (KIND=DP), INTENT(INOUT) :: state(2*varsize), xinit, xend
    REAL (KIND=DP), INTENT (IN) :: step
    INTEGER (KIND=DP), INTENT(IN) :: buffer

    REAL (KIND=DP) :: dstate(2*varsize)
    
    ! Temp variables to do evo: initial, end and current time in loop
    REAL (KIND=DP) :: tempinit, tempend, t_temp
    INTEGER (KIND=DP) :: i

    tempinit=xinit; tempend=xend

    CALL write_output(tempinit, 'time')

    ! Initial state may change slowly - run for some time without 
    ! checking for convergence:
    DO i=1, buffer
        t_temp = tempinit + step
        CALL evolve_state(state, tempinit, t_temp)
        CALL storage(BackCheck)
        tempinit = t_temp
        CALL write_output(tempinit, 'time')
     END DO

    ! Now check for convergence, exiting if reached steady state or max time:    
    time_evo_loop: DO WHILE(.NOT. converged())
       t_temp = tempinit + step
       CALL evolve_state(state, tempinit, t_temp)
       CALL storage(BackCheck)
       tempinit = t_temp
       CALL write_output(tempinit, 'time')
       
       ! If it hits maximum time limit:
       IF (tempinit .GE. tempend) THEN
          WRITE(51,'("# Warning: did not reach a steady state")')
          WRITE(52,'("# Warning: non-converged point")')
          EXIT time_evo_loop
       END IF

    END DO time_evo_loop
    
    ! Print out convergence info:
    IF (dbg_convergence) then
       DO i = 1, 10
          WRITE(*,*) StateStore(i, 1), StateStore(i, 2)
       END DO
    END IF
    
  END SUBROUTINE next_steady

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: evolve_state
  ! VARIABLES:
  !           state  - encoded state (INOUT)
  !           xinit  - initial state to pass to NAG (INOUT)
  !           xend   - final state from NAG (INOUT)
  !           eom    - subroutine to use as eqn of motion (INTERFACE)
  ! SYNOPSIS:
  !   This subroutine evolves the state from xinit to xend according 
  !   to eqns of motion.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE evolve_state(state, xinit, xend)
    USE nag_library, ONLY: d02bjf, d02bjw, d02bjx
    
    REAL (KIND=DP), INTENT(INOUT) :: state(2*varsize), xinit, xend

    REAL (KIND=DP) :: tol, W(40*varsize)
    INTEGER (KIND=DP) :: i
    INTEGER :: NN
    INTEGER :: ifail = 0

    ! NAG parameters:
    NN = 2*varsize; tol = 0.000001;

    !Now evolve the equation to the end point.
    !Output and stopping both suppressed.
    CALL d02bjf(xinit, xend, NN, state, diffeqn, tol, 'D', &
              & d02bjx, d02bjw, W, ifail)

  END SUBROUTINE evolve_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: storage
  ! VARIABLES:
  !           N - the number of elements stored (IN)
  ! SYNOPSIS:
  !   This subroutine stores the last N values of a chosen element 
  !   of the state. 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE storage(N)
    COMPLEX (KIND=DP) :: Phi(-nmax:nmax, -mmax:mmax), Lam
    INTEGER :: i, N
    
    CALL decode_state(CurrentState, Phi, Lam)
    
    DO i = 1, (N - 1)
       StateStore(i,:) = StateStore(i+1,:)
    END DO
    
    StateStore(N,1) = Real(Lam)
    StateStore(N,2) = Aimag(Lam)
    
  END SUBROUTINE storage

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: converged
  ! VARIABLES:
  ! SYNOPSIS:
  !   This function is true if the state has converged to steady state 
  !   (the last N derivatives are below tolerance).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION converged()
    LOGICAL :: converged
    INTEGER :: i
    REAL (KIND=DP) :: maximum_R, maximum_I, minimum_R, minimum_I
    
    converged = .TRUE.
    maximum_R = maxval(StateStore(:,1)); maximum_I = maxval(StateStore(:,2))
    minimum_R = minval(StateStore(:,1)); minimum_I = minval(StateStore(:,2))
    
    IF ( ( (maximum_R - minimum_R) .GT. convergence_threshold) .OR. &
         &( (maximum_I - minimum_I) .GT. convergence_threshold)) THEN
       converged = .FALSE.
    END IF
    
  END FUNCTION converged

END MODULE evolution
