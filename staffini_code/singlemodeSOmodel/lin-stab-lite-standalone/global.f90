MODULE global
  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER :: DP=KIND(1.0d0)
  
  COMPLEX (KIND=DP), PARAMETER :: i=(0.0D0,1.0D0)
  
  ! Physical parameters of the problem:
  REAL (KIND=DP) :: Kappa, N_A, Delta_a, E_0, omega, omega_r

  INTEGER (KIND=DP) :: lmax, imposed_nmax

CONTAINS 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: 
  !
  ! VARIABLES:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   interpolate
  ! VARIABLES:
  !           min, max - Range of variable
  !           i        - Index
  !           steps    - Number of points
  !           log      - Whether to use a log scale
  ! SYNOPSIS:
  !   Return a number interpolated between min and max, according to
  !   i/steps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION interpolate(min,max,i,steps,logarithmic)
    REAL (KIND=DP), INTENT(IN) ::min,max
    INTEGER, INTENT(IN) ::i,steps
    REAL (KIND=DP) ::interpolate
    LOGICAL :: logarithmic

    IF (steps .EQ. 1) THEN
       interpolate=min
    ELSE
       IF (logarithmic) THEN
          interpolate = min * exp( log(max/min)*(i-1)/(steps-1))
       ELSE
          interpolate=min + (max-min)*(i-1)/(steps-1)
       end IF
    end IF

    IF (i .EQ. 1) THEN
       interpolate = min
    ELSE IF (i .EQ. steps) THEN
       interpolate = max
    end IF
    

  END FUNCTION interpolate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   assert
  ! VARIABLES:
  !           condition
  !           errormsg
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE assert(condition,errormsg)
    LOGICAL, INTENT(IN) ::condition
    CHARACTER(LEN=*), INTENT(IN) ::errormsg

    IF (.NOT. condition) THEN
       WRITE(*,*) errormsg
       STOP
    end IF

  END SUBROUTINE assert

END MODULE global  
