MODULE sweeps
  USE global, ONLY: A_OmegaTil, A_Pump, Pump, OmegaTil, CurrentState, &
                    & tinit, tend, tstep, write_output, TimeOutPort, &
                    & TimeOutput, make_file
  USE evolution, ONLY: next_steady, diffeqn
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  ! This module contains all the subroutines necessary to perform sweeps.

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: sweep
  ! VARIABLES:
  !           var  - the variable in which direction to sweep (IN)
  !           dir  - sweep 'up' or 'down' in var direction (IN)
  !           min  - lower bound of sweep range (IN)
  !           max  - upper bound of sweep range (IN)
  !           no - number of steps in the sweep (IN) 
  ! SYNOPSIS:
  !   This subroutine sweeps up/down for either P or OmTil for a given 
  !   range by a given step.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE do_sweep(var, dir, min, max, no)
    REAL (KIND=DP), INTENT(INOUT) :: var
    REAL (KIND=DP), INTENT(IN) :: min, max, no
    CHARACTER(LEN=*) :: dir
    INTEGER (KIND=DP) :: buffer=1000

    SELECT CASE(dir)
    CASE ('up')
       var=min
       DO WHILE (var .LE. max) 
          WRITE(*,*) Pump, OmegaTil 
         
          ! Make a new time evo file to record evo:
          CALL make_file('time', TimeOutPort, TimeOutput)
          OPEN (TimeOutPort, file=TimeOutput, status='old', access='append')
          
          CALL next_steady(CurrentState,tinit, tend, tstep, buffer)
          CALL write_output(0.0d0, 'sweep')
          CLOSE(TimeOutPort)

          var=var+(max-min)/no
       END DO
    CASE ('down')
       var=max
       DO WHILE (var .GE. min)
          WRITE(*,*) Pump, OmegaTil

          ! Make a new time evo file to record evo:
          CALL make_file('time', TimeOutPort, TimeOutput)
          OPEN (TimeOutPort, file=TimeOutput, status='old', access='append')

          CALL next_steady(CurrentState, tinit,tend, tstep, buffer)
          CALL write_output(0.0d0, 'sweep')
          CLOSE(TimeOutPort)

          var=var-(max-min)/no
       END DO
    END SELECT

  END SUBROUTINE do_sweep

END MODULE sweeps
