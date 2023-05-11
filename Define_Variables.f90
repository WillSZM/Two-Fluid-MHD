module define_variables
    implicit none

    ! constants and normalizing parameters
    real(kind=8), parameter :: pi=3.1415926535
    real(kind=8), parameter :: q0=1.602176565d-19
    real(kind=8), parameter :: m0=9.10938356d-31
    real(kind=8), parameter :: eps0=8.854187817d-12
    real(kind=8), parameter :: mu0=1.2566370614d-6
    real(kind=8) :: x0
    real(kind=8) :: t0
    real(kind=8) :: n0
    real(kind=8) :: v0, E0, B0, pr0, Tem0, j0, R0, eta0
    real(kind=8) :: qe, qi
    real(kind=8) :: me, mi
    real(kind=8) :: const1, const2, const3

    
    ! grid and index
    integer :: Nx, Ny, Nz
    real(kind=8) :: xmax, xmin, ymax, ymin, zmax, zmin
    real(kind=8) :: hx, hy, hz
    real(kind=8), allocatable :: x(:), y(:), z(:)

    ! time and time step
    integer :: nstep, nmax, nout
    real(kind=8) :: time, tau
    real(kind=8), allocatable :: time_series(:)

    ! physical parameter
    real(kind=8), allocatable :: ne(:,:,:), ni(:,:,:)
    real(kind=8), allocatable :: vex(:,:,:), vey(:,:,:), vez(:,:,:), vix(:,:,:), viy(:,:,:), viz(:,:,:)
    real(kind=8), allocatable :: Bsx(:,:,:), Bsy(:,:,:), Bsz(:,:,:), Bex(:,:,:), Bey(:,:,:), Bez(:,:,:)
    real(kind=8), allocatable :: Esx(:,:,:), Esy(:,:,:), Esz(:,:,:), Eex(:,:,:), Eey(:,:,:), Eez(:,:,:)
    real(kind=8), allocatable :: Bsx_pre(:,:,:), Bsy_pre(:,:,:), Bsz_pre(:,:,:)
    real(kind=8), allocatable :: Esx_pre(:,:,:), Esy_pre(:,:,:), Esz_pre(:,:,:)
    real(kind=8), allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:), Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
    real(kind=8), allocatable :: jx(:,:,:), jy(:,:,:), jz(:,:,:)
    real(kind=8), allocatable :: Te(:,:,:), Ti(:,:,:)
    real(kind=8), allocatable :: pre(:,:,:), pri(:,:,:), eta(:,:,:)
    real(kind=8), allocatable :: Rex(:,:,:), Rey(:,:,:), Rez(:,:,:), Rix(:,:,:), Riy(:,:,:), Riz(:,:,:)
    real(kind=8), allocatable :: divB(:,:,:), divB_time(:), divE(:,:,:), divE_time(:)
    real(kind=8) :: divB_max, divE_max, gamma

    ! controling parameter
    integer :: cont
    logical :: is_abnormal_resistance


contains
    

end module define_variables
