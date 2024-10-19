PROGRAM linear_stability
  USE global, ONLY: i, interpolate, lmax, & 
                    & Kappa, omega_r, E_0, imposed_nmax

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  
  REAL (KIND=DP) :: Pmin, Pmax, P
  COMPLEX (KIND=DP) :: nu, omega_crit1, omega_crit2
  INTEGER :: iP, NP
  INTEGER, PARAMETER :: output=25
  CHARACTER :: filename*128

  ! Zeroes only at nu=0 (analytically):
  nu=0.0d0
 
  ! Range of P:
  READ(*,*) Pmin, Pmax, nP

  ! Read in model parameters:
  READ(*,*) E_0, Kappa, omega_r
  
  ! Read in number of levels of atoms:
  READ(*,*) lmax, imposed_nmax

  ! Setting up the file to write data:
  READ(*,*) filename
  OPEN(unit=output, file=filename, status='replace')

  ! Writing headers:
  WRITE(output, '("#",3X,A1,8X,2(2(A11,X),3X))') "P", "Re(w1)", "Im(w1)", &
       & "Re(w2)", "Im(w2)"
        
  PRINT*, "Finding transition at nu=0 directly..."
        
  DO iP=1, nP
     P=interpolate(Pmin, Pmax, ip, nP, .FALSE.)
     CALL find_omega_crit(P, omega_crit1, omega_crit2)
     
     ! Only saves the data with real solutions:
     IF (IMAGPART(omega_crit1) .LE. 1d-10 .AND. &
          & IMAGPART(omega_crit2) .LE. 1d-10) THEN
        WRITE(output,'(F8.3,4X,2(2(F11.2,X),3X))') P, omega_crit1, omega_crit2
     ELSE
     END IF
  END DO

  PRINT*, 'Written to file: ', filename

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: find_omega_crit
  ! Calculates critical omega(s) for a given P at nu=0.
  ! VARIABLES: P - input (IN)
  !            omega_crit1, 2 - first and second roots for each P
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

  SUBROUTINE find_omega_crit(P, omega_crit1, omega_crit2)
    USE nag_library, ONLY: c02ajf
    ! USE nag_f77_c_chapter, ONLY: c02ajf                             
    USE overlaps, ONLY:overlap_sum

    REAL (KIND=DP), INTENT(IN) :: P
    COMPLEX (KIND=DP), INTENT(OUT) :: omega_crit1, omega_crit2

    COMPLEX (KIND=DP) :: Sum_Term
    REAL (KIND=DP) :: A, B, C, root1(2), root2(2)
    INTEGER :: ifail

    ! Calling overlaps to calculate the terms of the l sum                                                                   
    CALL overlap_sum((0.0d0,0.0d0), P, Sum_Term) 
    
    ifail=1 

    A= -1
    IF (IMAGPART(Sum_Term) .EQ. 0) THEN
       B= 8 * E_0 * omega_r**2 * P**2 * REALPART(Sum_Term) 
    ELSE 
       PRINT*, "Sum term has an imaginary part"
    END IF

    C= -Kappa**2

    CALL c02ajf(A,B,C,root1,root2,ifail)
    
    omega_crit1=COMPLEX(root1(1), root1(2))
    omega_crit2=COMPLEX(root2(1), root2(2))
    
  END SUBROUTINE find_omega_crit

END PROGRAM linear_stability
