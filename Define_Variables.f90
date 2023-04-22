module define_variables
    !Define
    integer, parameter :: mx=101, my=101, mz=101

    logical :: lrstrt, uniformx, uniformy, uniformz, periody, cbndry, halfx, halfy, halfz
    logical :: is_foreta

    !about time and time step
    integer :: nstp(250), dnstep, nstop, nstep
    integer :: nend, nsp, ncase, nst, nst1, nint, npt, nst2
    real(kind=8) :: dtime, time, t0, dt

    !about grid and index
    integer :: nx, nxp1, ny, nyp1, nz2, nz2p1, nz, nzp1
    integer :: nsmthx, nsmthy, nsmthz
    real(kind=8) :: xmax, xmin, ymax, ymin, zmax, zmin, dxmin, dymin, dzmin
    real(kind=8) :: dx, dy, dz, ds
    real(kind=8) :: x(mx, my, mz, 8), xi(mx, my, mz, 8)
    real(kind=8) :: xm(mx, my, mz, 8)
    real(kind=8) :: xx(mx), yy(my), zz(mz)
    integer :: nnxh, nnx, nnxq, nnyh, nny, nnyq, nnzh, nnz, nnzq

    !physical parameter
    real(kind=8) :: gamma, di
    real(kind=8) :: vyi0
    real(kind=8) :: w0(mx, my, mz, 3), fx(mx, 8), fxi(mx, 8)
    real(kind=8) :: pr(mx, my, mz), etaf(mx, my, mz)
    real(kind=8) :: fs(mx, my, mz), gs(mx, my, mz), hs(mx, my, mz)
    real(kind=8) :: MagStr(mx, my, mz), j_para(mx, my, mz)

    integer :: cont
    real(kind=8) :: aw

    real(kind=8), parameter :: pi=3.1415926535
    real(kind=8) :: presc, hb2, hv2, pri

contains
    !Set Value
    subroutine set_values
        lrstrt = .false.
        uniformx = .true.
        uniformy = .true.
        uniformz = .false.
        periody = .false.
        halfx = .false.
        halfy = .false.
        halfz = .false.
        cbndry = .true.
        is_foreta = .false.

        ncase = 1
        nend = 7   !输出文件数量
        nst = 1
        nst1 = 1
        nst2 = 5
        cont = 2
        nint = 1

        dtime = 40.0    !每隔dtime输出一次
        dnstep = 1
        nstop = 100000

        nx = mx - 1
        ny = my - 1
        nz = mz - 1

        xmin = -15.d0
        xmax = 15.d0
        ymin = -50.d0
        ymax = 50.d0
        zmin = -50.d0
        zmax = 50.d0
        dxmin = 0.05d0
        dymin = 0.5d0
        dzmin = 0.5d0

        nsmthx = 2*mx/3
        nsmthy = 3*my/4
        nsmthz = 2*mz/3
        if (halfz) nsmthz = mz - 1

    end subroutine

end module define_variables
