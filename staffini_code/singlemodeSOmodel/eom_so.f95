MODULE eom_so
  USE global, ONLY: ii, encode_state, decode_state, N_A, nmax, mmax, & 
                  & OmegaTil, Pump, Phi, Lam, CurrentState, Kappa, &
                  & Omega_r, E_0, varsize
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: diffeqn
  ! VARIABLES:
  !           t - current time (IN)
  !           s - current (encoded) state (IN)
  !           ds_dt - (encoded) time derivative (OUT)
  ! SYNOPSIS:
  !   This takes current time t and state s and returns encoded ds_dt 
  !   according to the equations of motion derived in the notes.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE diffeqn(t, s, ds_dt)
    REAL (KIND=DP), INTENT(IN) :: t, s(*)
    REAL (KIND=DP), INTENT (OUT) :: ds_dt(*)
    REAL (KIND=DP) :: dummy

    COMPLEX (KIND=DP) :: dPhi(-nmax:nmax, -mmax:mmax), dLam
    
    INTEGER (KIND=DP) :: i,j,k,l

    dummy = t

    ! Decode encoded state for manipulation:
    CALL decode_state(s, Phi, Lam)
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Firstly, calculate the rate of change of Lam:
    
    dLam = (OmegaTil - ii*(Kappa/E_0)) * Lam
    
    ! Spanning all the cavity direction momenta (n)
    DO i = -nmax, nmax 
       
       ! Spanning all the pump direction momenta (m)
       DO j = -mmax, mmax
        
          ! Phi_(n+2s),m scattering term (k is s here):
          DO k = -1, 1, 2
             IF (ABS(i+2*k) .GT. nmax) THEN
                ! do nothing, ie. add zero to dLam
             ELSE
                dLam = dLam - Lam*CONJG(Phi(i+2*k,j))*Phi(i,j) 
             END IF
          END DO

          ! Phi_(n+s),(m+s') scattering term (k=s,l=s'):
          DO k = -1, 1, 2
             DO l = -1, 1, 2
                IF (ABS(i+k) .GT. nmax .OR. ABS(j+l) .GT. mmax) THEN
                   !do nothing, ie. add zero to dLam
                ELSE
                   dLam = dLam - Pump*CONJG(Phi(i+k,j+l))*Phi(i,j)
                END IF
             END DO
          END DO
    
       END DO

    END DO
    
    ! Add prefactors:
    dLam = -ii*E_0*dLam

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now calculate d_dt of the momentum matrix elements phi(i,j)
    
    ! Spanning cavity direction:
    DO i = -nmax, nmax
       
       ! Spanning entire pump direction:
       DO j = -mmax, mmax
          
          dPhi(i,j) = (i**2+j**2)*Phi(i,j)
          
        ! Phi_n,(m+2s) scattering term (k is s here):
          DO k = -1, 1, 2
             IF(ABS(j+2*k) .GT. mmax) THEN
                !do nothing
             ELSE
                dPhi(i,j) = dPhi(i,j) - (Pump**2)*Phi(i,j+2*k)
             END IF
          END DO

          ! Phi_(n+2s),m scattering term (k=s):                
          DO k = -1, 1, 2
             IF(ABS(i+2*k) .GT. nmax) THEN
                !do nothing
             ELSE
                dPhi(i,j) = dPhi(i,j) - (abs(Lam)**2)*Phi(i+2*k,j)
             END IF
          END DO
          
          ! Phi_(n+s),(m+s') scattering term (k=s,l=s'): 
          DO k = -1, 1, 2
             DO l = -1, 1, 2
                IF(ABS(i+k) .GT. nmax .OR. ABS(j+l) .GT. mmax) THEN
                   !do nothing
                ELSE
                   dPhi(i,j) = dPhi(i,j) - Pump*(Lam+conjg(Lam))*Phi(i+k,j+l)
                END IF
             END DO
          END DO
                
          ! Add prefactors:
          dPhi(i,j) = -ii*Omega_r*dPhi(i,j)
        
       END DO
    END DO
    
    !Encode this into the output state.
    call encode_state(dPhi, dLam, ds_dt)

  END SUBROUTINE diffeqn

END MODULE eom_so
