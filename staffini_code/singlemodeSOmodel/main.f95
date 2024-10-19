! It outputs into several files: the sweep file, which contains the
! steady state solutions at varying pump strength, with each run's 
! initial conditions being the steady state from the previous run, 
! and the time evolution of each indulvidual run also.

PROGRAM self_organisation
  USE global, ONLY: ii, Phi, Lam, CurrentState, nmax, mmax, varsize, Omega_r, &
       & A_Pump, A_OmegaTil, tinit, tend, tstep, E_0, Kappa, &
       & convergence_threshold, RunType, SweepVariable, N_A, make_file, &
       & SweepDir, Pump, OmegaTil, encode_state, imposed_nmax, lmax, &
       & TimeOutput, SweepOutput, TimeOutPort, SweepOutPort, decode_state

  USE evolution, ONLY: execute_single
  USE sweeps, ONLY: do_sweep
  USE linearstability, ONLY: OmegaTilCrit

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL (KIND=DP) :: norm, dPump, Cutoff, Saved_Pump_Start, counter, resolution, step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reading in parameters, setting up arrays:

  OPEN(12, file='params', status='old', action='read')

  READ(12,*) RunType, SweepVariable, SweepDir
  READ(12,*) nmax, E_0, Kappa, N_A, Omega_r
  READ(12,*) A_Pump, A_OmegaTil
  READ(12,*) tinit, tend, tstep
  READ(12,*) convergence_threshold
  READ(12,*) lmax, imposed_nmax

  CLOSE(12)

  mmax = nmax
  varsize = (2*nmax+1)*(2*mmax+1)+1

  ALLOCATE(CurrentState(2*varsize))
  ALLOCATE(Phi(-nmax:nmax, -mmax:mmax))

  ! Passing initial pump and omegatil
  Pump = A_Pump(1); OmegaTil = A_OmegaTil(1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up the initial state of the system:

  Phi = (1-ii)/(2*Sqrt(2*N_A)); Phi(0,0) = 1.0
  Phi(1,0)=0.0; Phi(0,1)=0.0; Phi(-1,0)=0.0; Phi(0,-1)=0.0
  Lam = 0.0
  norm = SUM(ABS(Phi)**2)
  Phi=Phi/sqrt(norm)

  CALL encode_state(Phi,Lam,CurrentState)
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Run program for specified runtype:

  SELECT CASE(RunType)

  CASE('omcrit')
     PRINT*, OmegaTilCrit(Pump)

  !~~~~~~~~~~~ TRICRITICAL POINT FINDING ~~~~~~~~~~~!

  CASE('tricrit')
     ! Takes a single pump value in and no omega

     CALL make_file('sweep', SweepOutPort, SweepOutput)
     OPEN(SweepOutPort, file=SweepOutput, status='old', access='append') 

     PRINT*, "Searching for tricritical point..."
     PRINT*, "At omega crit, sweeping P -"

     ! Set local omega as omega crit
     OmegaTil=OmegaTilCrit(Pump)

     ! Set dPump (interval on each side of stability line)
     dPump=0.1 ! originally 0.01

     ! Redefine pump intervals for sweeps:
     A_Pump(1)=Pump-dPump
     A_Pump(2)=Pump+dPump

     Saved_Pump_Start=A_Pump(1)
     Cutoff=0.1

     ! This keeps going right until lambda squared is suitably above noise:
     
     ! Standard sweep up:
     resolution=4.0d0
     CALL do_sweep(Pump, 'up', A_Pump(1), &
          & A_Pump(2), resolution)
     counter=resolution ! Keeping track of how many steps we're taking on the way up to come back down later
     step=2*dPump/resolution

     DO WHILE ( Abs(Lam)**2 .LT. Cutoff)
        ! Keep sweeping up bit by bit:
        A_Pump(1)=A_Pump(2)+step
        A_Pump(2)=A_Pump(2)+2*step
        CALL do_sweep(Pump, 'up', A_Pump(1), &
             & A_Pump(2), 1.0d0)
        counter=counter+2.0d0
        CALL decode_state(CurrentState, Phi, Lam)
     END DO
  
     ! Reset the start point of the pump:
     A_Pump(1)=Saved_Pump_Start

     ! Coming back down:
     CALL do_sweep(Pump, 'down', A_Pump(1), &
          & A_Pump(2), counter)

     PRINT*, "Sweep results written to:", SweepOutput

     CLOSE(SweepOutPort)

  !~~~~~~~~~~~~~~ SWEEP ROUTINES ~~~~~~~~~~~~~!   
  CASE('sweep')

     ! Make sweep file:
     CALL make_file('sweep', SweepOutPort, SweepOutput)
     OPEN(SweepOutPort, file=SweepOutput, status='old', access='append') 

     SELECT CASE(SweepVariable)

        ! Sweep omega tilde:
     CASE('omegatil')

        PRINT*, "Sweeping omegatil in direction: ", SweepDir
        CALL do_sweep(OmegaTil, SweepDir, A_OmegaTil(1), &
             & A_OmegaTil(2), A_OmegaTil(3))
        PRINT*, "Written to:", SweepOutput
        PRINT*, "Time evo written to:", TimeOutput

        ! Sweep P:   
     CASE('pump')

        PRINT*, "Sweeping P in direction: ", SweepDir
        CALL do_sweep(Pump, SweepDir, A_Pump(1), &
             & A_Pump(2), A_Pump(3))
        PRINT*, "Written to:", SweepOutput
        PRINT*, "Time evo written to:", TimeOutput
     END SELECT
     
     CLOSE(SweepOutPort)

  !~~~~~~~~~~~~~~ SINGLE-POINT TIME EVO ~~~~~~~~~~~~~!   
  CASE('single')

     ! Make time file only:
     CALL make_file('time', TimeOutPort, TimeOutput)
     OPEN(TimeOutPort, file=TimeOutput, status='old', access='append') 

     PRINT*, "Running time evolution for a single point in the phase diagram:"
     CALL execute_single()
     PRINT*, "Time evo written to:", TimeOutput

     CLOSE(TimeOutPort)

  END SELECT

  DEALLOCATE(Phi, CurrentState)

END PROGRAM self_organisation
