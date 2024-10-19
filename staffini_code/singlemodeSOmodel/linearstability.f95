MODULE linearstability
  USE global, ONLY: ii, Kappa, Omega_r, E_0
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0d0)

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
  ! NAME: OmegaTilCrit                                                                    
  ! VARIABLES:                                                                            
  !           OmegaTilCrit - critical omega til from linear stability                          
  !           P - input pump value (IN)                                             
  ! SYNOPSIS:                                                                             
  !   Returns the value of OmegaTilCrit for a given pump value.                            
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  FUNCTION OmegaTilCrit(P)
    REAL (KIND=DP) :: OmegaTilCrit
    REAL (KIND=DP), INTENT(IN) :: P
    
    COMPLEX (KIND=DP) :: OmegaSM, OmegaLG

    CALL find_omega_crit(P, OmegaSM, OmegaLG)

    IF (IMAGPART(OmegaLG) .LT. 1.0d-10) THEN
       OmegaTilCrit=REALPART(OmegaLG)/E_0
    ELSE
       OmegaTilCrit=0.0d0
       PRINT*, "Error: OmegaTilCrit is imaginary for the pump value."
    END IF

  END FUNCTION OmegaTilCrit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
  ! NAME: find_omega_crit  
  ! SYNOPSIS:
  !   Calculates critical omega(s) for a given P at nu=0.                                   
  ! VARIABLES: 
  !           P - input pump value (IN)                                                             
  !           omega_crit1, 2 - 1st and 2nd root for each P (OUT)                         
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

    ! Calling overlaps to calculate the terms of the l sum                               \
                                                                                          
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
  
END MODULE linearstability
