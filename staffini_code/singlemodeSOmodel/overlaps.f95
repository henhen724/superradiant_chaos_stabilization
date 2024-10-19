MODULE overlaps
  USE global, ONLY: ii

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=KIND(1.0D0)

  INTEGER (KIND=DP) :: iQ, NQ
  REAL (KIND=DP) :: overlap
  REAL (KIND=DP), ALLOCATABLE :: diag(:), offdiag(:)
  REAL (KIND=DP), ALLOCATABLE :: eigenvals(:), eigenvectors(:,:)
  REAL (KIND=DP), ALLOCATABLE :: gs_eigenvalue, gs_eigenvector(:) 
  REAL (KIND=DP), ALLOCATABLE :: es_eigenvector(:)

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: overlap_sum
  !
  ! Calculates the l sum term in the Det equation putting together 
  ! the eigenvalues and eigenvectors calculated by NAG routine
  !
  ! VARIABLES: nu - variable from main (IN)
  !            P - variable from main (IN)
  !            Sum_Term - value of computed l sum for Det (OUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE overlap_sum(nu, P, Sum_Term)
    USE global, ONLY: Omega_r, lmax

    REAL (KIND=DP), INTENT(IN) :: P
    COMPLEX (KIND=DP), INTENT(IN) :: nu
    COMPLEX (KIND=DP), INTENT(OUT) :: Sum_Term
    REAL (KIND=DP) :: Delta_epsilon

    ! EVEN SECTOR
    CALL generate_h(.FALSE.,P,diag,offdiag)
    CALL get_eigenstuff(diag, offdiag,eigenvals,eigenvectors)
    
    ALLOCATE(gs_eigenvector(SIZE(eigenvectors,1)))
    gs_eigenvector=eigenvectors(:,1)
    gs_eigenvalue=eigenvals(1)
    
    DEALLOCATE (diag, offdiag, eigenvals, eigenvectors)
    
    !ODD SECTOR
    CALL generate_h(.TRUE.,P,diag,offdiag)
    CALL get_eigenstuff(diag, offdiag,eigenvals,eigenvectors)
    
    ALLOCATE(es_eigenvector(SIZE(eigenvectors,1)))
    
    NQ=MIN(SIZE(eigenvectors, 2),lmax)
    Sum_Term=0
    
    DO iQ=1,NQ
       
       es_eigenvector=eigenvectors(:,iQ)
       
       CALL get_overlap(gs_eigenvector,es_eigenvector,overlap)
       
       ! WRITE(out_unit,'(X,D12.5,X,I12,X,3(D12.5,X))') P, iQ, &
       !      & gs_eigenvalue, eigenvals(iq), overlap
       
       ! Define for convenience
       Delta_epsilon=-gs_eigenvalue+eigenvals(iq)+1.0D0
       
       ! Equation of the sum term from notes (double check)
       Sum_Term=Sum_Term + &
            & (Delta_epsilon*overlap**2)/((omega_r*Delta_epsilon)**2-nu**2)
       
    END DO

    DEALLOCATE (diag, offdiag, eigenvals,eigenvectors)
    DEALLOCATE(gs_eigenvector, es_eigenvector)
    
  END SUBROUTINE overlap_sum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE: generate_h
  !
  ! Builds diagonal and off-diagonal arrays for the  tridiagonal 
  ! hamiltonian matrix for either even or odd subspace.
  !
  ! VARIABLES: is_odd - logical toggle (IN)
  !            P - from main, matrix size is then 100P(+1) (IN)
  !            diag - the diagonal entries of h (OUT)
  !            offdiag - the offdiagonal entries of h (OUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE generate_h(is_odd,P,diag,offdiag)
    USE global, ONLY: imposed_nmax
    REAL (KIND=DP) :: P
    LOGICAL, INTENT(IN) :: is_odd
    REAL (KIND=DP), INTENT(OUT), ALLOCATABLE :: diag(:), offdiag(:)
    INTEGER (KIND=DP) :: XX,A, max_n, nmax
 
    ! A is the size of the large matrix we want to use in multilevel:
    A=50

    ! then imposed_nmax allows to truncate the matrix to a true two-level
    ! Dicke treatment.

    IF (is_odd) THEN
       ! For odd matrix, want either 2x2 (Dicke) or odd-sized:
       nmax=2*INT(A*P)+11
       max_n=MIN(imposed_nmax, nmax)
       ALLOCATE (diag(max_n+1))
       diag=(/ (XX**2, XX=-max_n,max_n, 2) /)
    ELSE 
       nmax=2*INT(A*P)+10
       ! For even, want either single state (1x1) or even-sized:
       max_n=MIN(imposed_nmax-1, nmax)
       ALLOCATE (diag(max_n+1))
       diag=(/ (XX**2, XX=-max_n,max_n, 2) /)
    END IF
    
    ALLOCATE (offdiag(max_n))
    offdiag=-P**2
    
  END SUBROUTINE generate_h

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SUBROUTINE: get_eigenstuff
  !
  ! Calls NAG routines to compute eigenvalues and eigenvectors. 
  !
  ! VARIABLES: dg - diagonal of the tridiagonal matrix (IN)
  !            od - off-diagonal of tridiagonal matrix (IN)
  !            evals - computed eigenvalues (in an array) (OUT)
  !            evects - computed eigenvectors (in 2d array) (OUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_eigenstuff(dg,od,evals,evects)
    USE nag_library, ONLY: f08jjf, f08jkf
    !USE nag_f77_f_chapter, ONLY: f08jjf, f08jkf

    REAL (KIND=DP), INTENT(IN) ::  dg(:), od(:)
    REAL (KIND=DP), ALLOCATABLE, INTENT(OUT) :: evals(:), evects(:,:)

    CHARACTER(1), PARAMETER :: range='A', order='B'
    REAL (KIND=DP), PARAMETER :: abstol=0.0

    INTEGER :: il, iu, m, nsplit, info, dim
    INTEGER, ALLOCATABLE :: iwork(:), iblock(:), isplit(:), ifailv(:)
    REAL (KIND=DP) :: vl, vu
    REAL (KIND=DP), ALLOCATABLE :: work(:)

    dim=SIZE(dg)

    ALLOCATE (iblock(dim), isplit(dim),iwork(3*dim), work(5*dim), ifailv(dim))
    ALLOCATE (evals(dim),evects(dim,dim))

    !   uses NAG routine to get eigenvals of h
    CALL f08jjf(range, order, dim, vl, vu, il, iu, abstol, dg, &
         & od, m, nsplit, evals, iblock, isplit, & 
         & work, iwork, info)

    !   uses NAG routine and results from get_eigenvals
    !   to get eigenvectors of h

    CALL f08jkf(dim, dg, od, dim, evals, iblock, isplit, evects, dim, work, &
         & iwork, ifailv, info)

    DEALLOCATE(iwork,work,ifailv,iblock,isplit)

  END SUBROUTINE get_eigenstuff

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME: get_overlap
  !
  ! Calculates the overlap term of the equation.
  !
  ! VARIABLES: ev_gs - eigenvector of ground state (IN)
  !            ev_es - eigenvector of excited state (IN)
  !            overlap - computed overlap (OUT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_overlap(ev_gs,ev_es,overlap)
    USE global, ONLY: assert

    REAL (KIND=DP), INTENT(IN) ::ev_gs(:),ev_es(:)
    REAL (KIND=DP), INTENT (OUT) :: overlap
    INTEGER :: dim

    dim=SIZE(ev_gs)
    CALL assert(SIZE(ev_es).EQ.dim+1, "Wrong size of excited state array")

    overlap= SUM(ev_gs(1:dim) * ev_es(1:dim)) &
         & + SUM(ev_gs(1:dim) * ev_es(2:dim+1)) 

  END SUBROUTINE get_overlap

END MODULE overlaps
