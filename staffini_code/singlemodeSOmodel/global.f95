MODULE global
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  COMPLEX (KIND=DP), PARAMETER :: ii=(0.0D0, 1.0D0)

  ! Global state variables to update:
  COMPLEX (KIND=DP), ALLOCATABLE :: Phi(:,:)
  REAL (KIND=DP), ALLOCATABLE :: CurrentState(:)
  COMPLEX (KIND=DP) :: Lam
  REAL (KIND=DP) :: OmegaTil, Pump

  ! Flags defining program operation
  CHARACTER (len=15) :: RunType, SweepVariable, SweepDir

  ! Physical parameters of the problem:
  REAL (KIND=DP) :: E_0, Kappa, N_A, Omega_r 

  !Arrays to hold range on parameters pump and omegatil (min, max and number)
  REAL (KIND=DP) :: A_Pump(3), A_OmegaTil(3)

  ! Initial time, max time cutoff and time step for each run
  REAL (KIND=DP) :: tinit, tend, tstep

  ! x and y size of momentum-space matrix and no. of complex vars in problem
  INTEGER :: nmax, mmax, varsize

  ! Convergence threshold:
  REAL (KIND=DP) :: convergence_threshold

  ! Output files:
  CHARACTER (len=100) :: TimeOutput, SweepOutput
  INTEGER (KIND=DP), PARAMETER :: TimeOutPort=51, SweepOutPort=52

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! For linear stability:
  INTEGER (KIND=DP) :: lmax, imposed_nmax

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: encode_state
  ! VARIABLES:
  !           Phi   - the phi momentum matrix (IN)
  !           Lam   - the light field (IN)
  !           state - the resulting state vector (OUT)
  ! SYNOPSIS:
  !   The routine flattens the matrix into a state vector for the
  !   NAG routines to handle.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE encode_state(Phi, Lam, state)
    COMPLEX (KIND=DP), INTENT(IN) :: Phi(:,:), Lam
    REAL (KIND=DP), INTENT(OUT) :: state(*)
    INTEGER (KIND=DP) :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    state(1:siz) = reshape(real(Phi), (/ siz /))
    state(siz+1:2*siz) = reshape(aimag(Phi),(/ siz /))
    state(2*siz+1:2*varsize) = (/ real(Lam), aimag(Lam) /)
    
  END SUBROUTINE encode_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: decode_state
  ! VARIABLES:
  !           state - the encoded state vector (IN)
  !           Phi   - the phi momentum matrix (OUT)
  !           Lam   - the light field (OUT)
  ! SYNOPSIS:
  !   Extracts phi matrix and scalar lambda from state vector.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE decode_state(state, Phi, Lam)
    COMPLEX (KIND=DP), INTENT(OUT) :: Phi(:,:), Lam
    REAL (KIND=DP), INTENT(IN) :: state(*)
    INTEGER (KIND=DP) :: siz
    
    siz = (2*nmax+1)*(2*mmax+1)
    
    Phi =    reshape(state(1:siz), (/ 2*nmax+1, 2*mmax+1 /)) + &
          & ii*reshape(state(siz+1:2*siz), (/ 2*nmax+1, 2*mmax+1 /))
          
    Lam = state(2*siz+1) + ii*state(2*varsize)
    
  END SUBROUTINE decode_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: write_output
  ! VARIABLES:
  !           t    - the physical time t (IN)
  !           type - the type of run (pump sweep or time evolution) (IN)
  ! SYNOPSIS:
  !   This subroutine writes the relevant information to file.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_output(t, type)
    REAL (KIND=DP), INTENT(IN) :: t
    CHARACTER (len=*), INTENT(IN) :: type
    
    real (KIND=kind(1.0D0)) :: sz, ipr
    complex (KIND=kind(1.0D0)) :: sp

    CALL decode_state(CurrentState, Phi, Lam)
    
    ! Dicke spin components:
    sz=0.5*(abs(2*Phi(1,1))**2 - abs(Phi(0,0))**2) ! S^z
    sp=2*Conjg(Phi(1,1))*Phi(0,0) ! S+

    ! Inverse participation ratio:
    ipr = SUM(ABS(phi)**4) / ( (SUM(ABS(phi)**2))**2 )
    
    ! Stuff to write out depending on type of run:
    SELECT CASE(type)
    CASE('time')
       WRITE(51,'(X,6(D14.7,X))') t/1000.0, Lam, ABS(phi(0,0))**2, &
            & ipr, ABS(Phi(1,1))**2
    CASE('sweep')
       WRITE(52,'(X,9(D14.7,X))') Pump, OmegaTil, Abs(Lam)**2, &
            & Lam, sp, sz, ipr
    END SELECT
        
  END SUBROUTINE write_output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: dump_params
  ! VARIABLES:
  !           handle - the output file to write to (IN)
  ! SYNOPSIS:
  !   Writes the parameters used to a file for good measure.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE dump_params(handle)
    integer (KIND=DP), intent(in) :: handle

    write(handle,'("# Size of grid:",I5,I5, "Num. Atoms",D12.5)') &
         & nmax, mmax, N_A
    write(handle,'("# E0=",D12.5," Kappa=",D12.5," OmegaR=",D12.5)') &
         & E_0, Kappa, Omega_r
    write(handle,'("# Initial OmegaTil=",D12.5," Pump=",D12.5)') &
         & OmegaTil, Pump

  END SUBROUTINE dump_params 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: make_file
  ! VARIABLES:
  !           
  ! SYNOPSIS:
  !   Sets up output files and writes headers depending on sweep type.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_file(type,port,namestring)
    CHARACTER(LEN=*), INTENT(IN) :: type
    INTEGER(KIND=DP), INTENT(IN) :: port
    CHARACTER(LEN=*), INTENT(OUT) :: namestring

    SELECT CASE(type)

  !~~~~~~~~~~ Create file for time evolution data:
    CASE ('time')

       ! Write the filename variable:
       WRITE(namestring,"(a9, f5.3, a10, f5.3, a4)") &
            & "Evo-Pump_", Pump, "-OmegaTil_", OmegaTil, ".dat"
       namestring=TRIM(namestring)

       ! Create the file:
       OPEN(port, file=namestring, status='replace')

       ! Dump parameters and write the headers:
       CALL dump_params(port)
       WRITE(port,'("#",6(A14,X))') "Time(ms)", "R[Lam]", "I[Lam]", &
               &   "|Phi(0,0)|^2", "IPR", "|Phi(1,1)|^2"

       ! Close the file:
       CLOSE(port)

  !~~~~~~~~~~ Create file for sweep data:
    CASE ('sweep')
       
       ! Write the filename variable:
       WRITE(namestring,"(a15, f5.3, a4, f5.3, a3, f5.3, a4, f5.3, a4)") &
            & "Sweep-OmegaTil_", A_OmegaTil(1), "_to_", A_OmegaTil(2), &
            & "-P_", A_Pump(1), "_to_", A_Pump(2), ".dat"
       namestring=TRIM(namestring)
       
       ! Create the file:
       OPEN(port, file=namestring, status='replace')

       ! Dump parameters and write the headers:
       CALL dump_params(port)
       WRITE(port,'("#",9(A14,X))') "Pump","OmegaTil", "|Lam|^2", "R[Lam]", "I[Lam]", &
             &   "R[Sp]", "I[Sp]", "Sz", "IPR"

       ! Close the file:
       CLOSE(port)

    END SELECT
    
  END SUBROUTINE make_file
      

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
  ! NAME:   assert                                                                        
  ! VARIABLES:                                                                            
  !           condition                                                                   
  !           errormsg                                                                    
  ! SYNOPSIS:                                                                             
  !   Returns errormsg if condition fails.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                     

  SUBROUTINE assert(condition,errormsg)
    LOGICAL, INTENT(IN) ::condition
    CHARACTER(LEN=*), INTENT(IN) ::errormsg

    IF (.NOT. condition) THEN
       WRITE(*,*) errormsg
       STOP
    end IF

  END SUBROUTINE assert

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

END MODULE global


