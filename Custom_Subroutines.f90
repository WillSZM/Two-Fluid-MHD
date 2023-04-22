module custom_subroutines
    use define_variables
    use custom_functions
    implicit none

contains

    subroutine input
        !
        ! --------------------------
        !  This routine inputs parameters and define basic variables
        !   LRSTRT: =.f., starting from t=0; =.t., continueing from
        !           steps as given by NST.
        !   NEND:   the final steps intended for the current run ,
        !           including the steps from the previous run.
        !   NSP:    time step interval for diagnostic plots.
        !   DELTAV: thickness of the initial velocity shear.
        !   DELTAC: thickness of the initial current sheet (for
        !           normalization coventions , see routine INITIA.
        !   RATIO:  ratio of dy to dx, or the ratio of the width to
        !           length of the box, if the number of grid points is such
        !           that my=mx.
        !   NSMTHX:  the number of rows starting from incoming
        !           boundary for the smoothing region in x-direction.
        !   NSMTHY:  the number of rows starting from incoming
        !           boundary for the smoothing region in y-direction.
        !   ETA:    exact inverse of magnetic Renolds number.
        !   GAMMA:  adiabatic constant.
        !   BETAS:   ratio of kinetic to magnetic pressure at magnetosheath side
        !   MACHS:   Alfven mach number of incoming flow  at magnetosheath side.
        !   MACHM:   Alfven mach number of incoming flow at magnetopause side.
        !   ANGS:    bz=b*cos(ang) and by=b*sin(ang), at magnetosheath side.
        !   ANGM:    bz=b*cos(ang) and by=b*sin(ang), at magnetopause side.
        !   TS0:   initia magnetosheath temperature
        !    tm0:   initia magnetopause temperature.
        !   BS0:   initia magnetosheath magnetic field strength
        !   BM0:   initia magnetopause magnetic field strength.
        !   NCASE:  case number of the run.
        !   NST:    beginning data file number for the current run.
        !   NINT:   data file number increment.
        ! --------------------------
        !

        ! Coordinates system
        !                A z     /
        !                l     /
        !                l   /
        !                l /
        !  <-------------*--------------
        !  x            /l
        !             /  l
        !           /    l
        !         /      l
        !        y
        implicit none

        integer :: jx, jz

        ! nx = mx - 1
        ! ny = my - 1
        ! nz = mz - 1

!         call gridpnt
!         open (unit=11, file='grid.dat', status='unknown', form='formatted')
!         write (11, 99) (xx(jx), jx=1, mx), (zz(jz), jz=1, mz)
! 99      format(5(1x, e10.4))

        ! time = 0.0
        ! nstep = 0
        ! nst1 = 1

        return
    end subroutine

    subroutine initia

        !----------------
        ! Defines coordinates system and specifies initial configuration.
        ! Normlization convention:
        !   1. Density --- normalised to asymtotic value, i.e., rho=1
        !   2. Magnetic field --- normalised to asymtotic value, i.e.,
        !                         b0=1.
        !   3. Velocity --- normalised to asymtotic Alfven speed, VA=1, a
        !                   natural result of 1. and 2.
        !   4. Length --- normalised to a=10*dx, i.e., dx=0.1
        !   5. Time --- normalised to a/VA.
        !---------------
    
        ! assign asymmetric quantities:

        implicit none

        real(kind=8) :: rhom, rhos, pp, gaminv, phi, phirad, delbz, by0, bz0, pmsp, pmsh, betas, p0, bzm, bzs
        real(kind=8) :: betam, bm0, bs0, epsilon
        real(kind=8) :: signy, thita_Xline, thita_Xline_rad, tan_Xline
        integer :: jx, jy, jz

        time = 0.0
        t0 = 0.0
        nstep = 0

        gamma = 1.66667d0
        di = 0.3d0
        aw = 1d0    !used in function foreta

        betam = 0.01d0  !beta in magnetopause
        bm0 = 1.d0  !initia magnetopause magnetic field strength
        bs0 = 1.d0  !initia magnetosheath magnetic field strength
        pmsp = betam*(0.5*bm0**2)   !thermal pressure in magnetopause
        pmsh = pmsp + 0.5*(bm0**2 - bs0**2) !thermal pressure in magnetosheath
        betas = pmsh/(0.5*bs0**2)   !beta in magnetosheath

        epsilon = 0.d0
        rhom = 0.1
        rhos = 0.1
        pp = 0.
        gaminv = 1./gamma

        phi = 180.d0
        phirad = 2.*pi*phi/360.
        delbz = 0.5*sqrt(bm0**2 + bs0**2 - 2.*bm0*bs0*cos(phirad))
        by0 = 0.5*bs0*sin(phirad)/delbz
        bz0 = 0.25*(bm0**2 - bs0**2)/delbz
        bzm = bz0 + delbz
        bzs = bz0 - delbz

        vyi0 = 0.d0

        signy = 1
        if (phi .lt. 180) signy = -1

        !cjx  the angle between the X-line(Y-direction) and MR plane
        thita_Xline = 0.0
        thita_Xline_rad = 2.*pi*thita_Xline/360.
        tan_Xline = tan(thita_Xline_rad)
        !cjx  -------------------------------------------------

        ! xi is exactly the same with x in initialization, so just make them equal
        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            x(jx, jy, jz, 1) = 1.0
            
            x(jx, jy, jz, 2) = 0.0
            x(jx, jy, jz, 3) = 0.0
            x(jx, jy, jz, 4) = 0.0

            x(jx, jy, jz, 5) = 0.0
            x(jx, jy, jz, 6) = 0.0
            x(jx, jy, jz, 7) = -tanh(xx(jx))

            etaf(jx, jy, jz) = 0.0d0

        end do
        end do
        end do
        !cjx----------------------------------------------------------------
        !      do 3 jz=1,mz
        !      do 3 jy=1,my
        !      do 3 jx=1,mx
        !      x(jx,jy,jz,1) = 0.5*(rhom+rhos) +0.5*(rhom-rhos)*tanh(xx(jx)/aw)
        !      x(jx,jy,jz,5)  = 0.0
        !      x(jx,jy,jz,6)  = signy*sqrt(by0*by0+epsilon*delbz*delbz
        !     1              /cosh(xx(jx)/aw)**2)
        !      x(jx,jy,jz,7)  = bz0 - delbz*tanh(xx(jx)/aw)
        !      xi(jx,jy,jz,1) = 0.5*(rhom+rhos)
        !     1   +0.5*(rhom-rhos)*tanh((xx(jx)+dx/2.)/aw)
        !      xi(jx,jy,jz,5)  = 0.0
        !      xi(jx,jy,jz,6)  = signy*sqrt(by0*by0+epsilon*delbz*delbz
        !     1              /cosh((xx(jx)+dx/2.)/aw)**2)
        !      xi(jx,jy,jz,7)  = bz0 - delbz*tanh((xx(jx)+dx/2.)/aw)
        !    3 continue

        !      do 4 jz=1,mz
        !      do 4 jy=1,my
        !      do 4 jx=1,mx
        !      x(jx,jy,jz,2)=0.
        !      x(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*(1.-tanh(xx(jx)/aw))
        !     1                    *x(jx,jy,jz,1)/bs0
        !      x(jx,jy,jz,4)=0.5*v0*x(jx,jy,jz,7)*(1.-tanh(xx(jx)/aw))
        !     1                    *x(jx,jy,jz,1)/bs0
        !      xi(jx,jy,jz,2)=0.
        !      xi(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*
        !     1        (1.-tanh((xx(jx)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
        !      xi(jx,jy,jz,4)  = 0.5*v0*x(jx,jy,jz,7)*
        !     1        (1.-tanh((xx(jx)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
        !    4 continue
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        presc = 1.0 !total pressure in magnetopause

        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            hb2 = .5*(x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2)
            hv2 = .5*(x(jx, jy, jz, 2)**2 + x(jx, jy, jz, 3)**2 + x(jx, jy, jz, 4)**2)/x(jx, jy, jz, 1)
            pr(jx, jy, jz) = presc - hb2
            x(jx, jy, jz, 8) = hv2 + hb2 + pr(jx, jy, jz)/(gamma - 1)
        end do
        end do
        end do

        xi = x
        
        call current(xi, 1)

        xi(:,:,:,3) = xi(:,:,:,3) + vyi0*di*w0(:,:,:,2)

        call current(x, 1)

        x(:,:,:,3) = x(:,:,:,3) + vyi0*di*w0(:,:,:,2)

        fx = x(:,1,1,:)
        fxi = xi(:,1,1,:)

        xm = 0
        do jz = 1, mz
            do jy = 1, my
                xi(:,jy,jz,:) = x(:,jy,jz,:) - fxi
            end do
        end do

        return
    end subroutine

    subroutine stepon

        !     This routine time-advances X's by 2 step Lax-Wendroff scheme
        !     note: X is always the up-to-date value while Xm being the
        !           intermediate value.


        implicit none

        real(kind=8) :: hdt, cmax ,cmin, xcmax, ycmax, zcmax, xcmin, ycmin, zcmin, caf0
        real(kind=8) :: fdx(mx,my,mz)=0, gdy(mx,my,mz)=0, hdz(mx,my,mz)=0, ppx(mx,my,mz)=0
        integer :: jx, jy, jz, m, k
        !cmax and cmin seems not been used
        !
        ! 1. first step
        !
        ! 1.1 Calculate fluxes

        ! 1.2 Advance the first step
        hdt = 0.5*dt
        call current(x, 1)
        call foreta(time, 1)
        call pressure(x, 1)

        do m = 1, 8

            call flux(x, fs, gs, hs, mx, my, mz, m, 1)

            fdx(1:nx,1:ny,1:nz) = -hdt*( fs(2:nx+1,2:ny+1,2:nz+1) - fs(1:nx,2:ny+1,2:nz+1) &
                                        + fs(2:nx+1,2:ny+1,1:nz) - fs(1:nx,2:ny+1,1:nz) &
                                        + fs(2:nx+1,1:ny,2:nz+1) - fs(1:nx,1:ny,2:nz+1) &
                                        + fs(2:nx+1,1:ny,1:nz) - fs(1:nx,1:ny,1:nz) )/(4.0*dx)
            gdy(1:nx,1:ny,1:nz) = -hdt*( gs(2:nx+1,2:ny+1,2:nz+1) - gs(2:nx+1,1:ny,2:nz+1) &
                                        + gs(2:nx+1,2:ny+1,1:nz) - gs(2:nx+1,1:ny,1:nz) &
                                        + gs(1:nx,2:ny+1,2:nz+1) - gs(1:nx,1:ny,2:nz+1) &
                                        + gs(1:nx,2:ny+1,1:nz) - gs(1:nx,1:ny,1:nz) ) / (4.0*dy)
            hdz(1:nx,1:ny,1:nz) = -hdt*( hs(2:nx+1,2:ny+1,2:nz+1) - hs(2:nx+1,2:ny+1,1:nz) &
                                        + hs(2:nx+1,1:ny,2:nz+1) - hs(2:nx+1,1:ny,1:nz) &
                                        + hs(1:nx,2:ny+1,2:nz+1) - hs(1:nx,2:ny+1,1:nz) &
                                        + hs(1:nx,1:ny,2:nz+1) - hs(1:nx,1:ny,1:nz) ) / (4.0*dz)
            ppx(1:nx,1:ny,1:nz) = ( xi(1:nx,1:ny,1:nz,m) + xi(2:nx+1,1:ny,1:nz,m) + xi(1:nx,2:ny+1,1:nz,m) &
                                    + xi(1:nx,1:ny,2:nz+1,m) + xi(2:nx+1,2:ny+1,1:nz,m) + xi(1:nx,2:ny+1,2:nz+1,m) &
                                    + xi(2:nx+1,1:ny,2:nz+1,m) + xi(2:nx+1,2:ny+1,2:nz+1,m) ) / 8.0

            xm(1:nx,1:ny,1:nz,m) = fdx(1:nx,1:ny,1:nz) + gdy(1:nx,1:ny,1:nz) &
                                    + hdz(1:nx,1:ny,1:nz) + ppx(1:nx,1:ny,1:nz)
            do jz = 1, nz
                do jy = 1, ny
                    xm(:,jy,jz,m) = xm(:,jy,jz,m) + fxi(:,m)
                end do
            end do

        end do

        
        
        
        ! 2. Second step
        
        ! 2.1 Calculate fluxes
        ! 2.2 Advance the second time step
        call current(xm, 2)
        call foreta(time, 2)
        call pressure(xm, 2)

        do m = 1, 8
            call flux(xm, fs, gs, hs, nx, ny, nz, m, 2)

            fdx(2:nx,2:ny,2:nz) = -dt*( fs(2:nx,2:ny,2:nz) - fs(1:nx-1,2:ny,2:nz) &
                                        + fs(2:nx,2:ny,1:nz-1) - fs(1:nx-1,2:ny,1:nz-1) &
                                        + fs(2:nx,1:ny-1,2:nz ) - fs(1:nx-1,1:ny-1,2:nz ) &
                                        + fs(2:nx,1:ny-1,1:nz-1) - fs(1:nx-1,1:ny-1,1:nz-1) ) / (4.0*dx)
            gdy(2:nx,2:ny,2:nz) = -dt*( gs(2:nx,2:ny,2:nz) - gs(2:nx ,1:ny-1,2:nz ) &
                                        + gs(2:nx,2:ny,1:nz-1) - gs(2:nx,1:ny-1,1:nz-1) &
                                        + gs(1:nx-1,2:ny,2:nz) - gs(1:nx-1,1:ny-1,2:nz) &
                                        + gs(1:nx-1,2:ny,1:nz-1) - gs(1:nx-1,1:ny-1,1:nz-1) ) / (4.0*dy)
            hdz(2:nx,2:ny,2:nz) = -dt*( hs(2:nx,2:ny,2:nz) - hs(2:nx,2:ny,1:nz-1) &
                                        + hs(2:nx,1:ny-1,2:nz) - hs(2:nx,1:ny-1,1:nz-1) &
                                        + hs(1:nx-1,2:ny,2:nz) - hs(1:nx-1,2:ny,1:nz-1) &
                                        + hs(1:nx-1,1:ny-1,2:nz) - hs(1:nx-1,1:ny-1,1:nz-1) ) / (4.0*dz)
            x(2:nx,2:ny,2:nz,m) = x(2:nx,2:ny,2:nz,m) + fdx(2:nx,2:ny,2:nz) + gdy(2:nx,2:ny,2:nz) + hdz(2:nx,2:ny,2:nz)

        end do

        call bndry(x, 1)

        !      if(mod(nstep,10).eq.0) call cleanb
        caf0 = 1.-0.25*tanh(time/100.0)

        do jz = 1, mz
        do jy = 1, my
            xi(:, jy, jz, :) = x(:, jy, jz, :) - fx(:, :)
        end do
        end do


        !      call smthf(xi,caf0)

        if (mod(nstep, 3) .eq. 0) then
            call smthxyz(xi, 0.875d0, 1)    !smooth for the first time
            call smthxyz(xi, 0.996d0, 1)    !smooth for the second time
        end if

        do jz = 1, mz
        do jy = 1, my
            x(:, jy, jz, :) = xi(:, jy, jz, :) + fx(:, :)
        end do
        end do


        if (mod(nstep, 10) .eq. 0) then
            call current(x, 1)
            do jz = 1, mz
            do jy = 1, my
            do jx = 1, mx
                ! magnetic field strength
                MagStr(jx, jy, jz) = sqrt(x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2)
                ! j parallel
                j_para(jx, jy, jz) = (w0(jx, jy, jz, 1)*x(jx, jy, jz, 5) + w0(jx, jy, jz, 2) &
                                     *x(jx, jy, jz, 6) + w0(jx, jy, jz, 3)*x(jx, jy, jz, 7))/MagStr(jx, jy, jz)
            end do
            end do
            end do

            open (unit=14, file='facmax.dat', status='unknown', form='formatted')
            write (14, 12) time, nstep, ncase
            k = 0
            do while (k < 8)
                k = k + 1
                jz = (mz*k/9.) + 1
                cmax = -10.
                cmin = 10.

                do jy = 1, my
                do jx = 1, mx - 2
                    if (xm(jx, jy, jz, 2) .gt. cmax) then
                        cmax = j_para(jx, jy, jz)
                        xcmax = xx(jx)
                        ycmax = yy(jy)
                        zcmax = zz(jz)
                    end if
                    if (xm(jx, jy, jz, 2) .lt. cmin) then
                        cmin = j_para(jx, jy, jz)
                        xcmin = xx(jx)
                        ycmin = yy(jy)
                        zcmin = zz(jz)
                    end if
                end do
                end do

                write (14, 13) xcmin, ycmin, zcmin, xcmax, ycmax, zcmax, cmin, cmax
            end do

12          format(1x, f8.2, 2(2x, i5))
13          format(6(1x, f6.1), 2(1x, e11.3))
            !      if(mod(nstep,30).eq.0) then
            !      open(unit=15,file='sat.dat',status='unknown',form='formatted')
            !      write(15,12)(time,nstep,ncase)
            !      write(15,17)((xm(jx,jy,6,2),jx=11,nx,10),jy=1,my,4)
            !   17 format(8(1x,e9.3))
            !      endif
        end if

        if (mod(nstep, 20) .eq. 0) call energy

        return
    end subroutine

    subroutine facur
        implicit none

        character*8 output

        real(kind=8) :: fact(mx, my, mz), facp(mx, my, mz), faci(mx, my, mz)
        real(kind=8) :: pp1, pp2, pp3, pp4, pp5, pp6
        integer :: jx, jy, jz, m, k

        call current(x, 1)
        call pressure(x, 1)

        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            xm(jx, jy, jz, 1) = sqrt(x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2)
            xm(jx, jy, jz, 2) = (w0(jx, jy, jz, 1)*x(jx, jy, jz, 5) + w0(jx, jy, jz, 2)*x(jx, jy, jz, 6) &
                                 + w0(jx, jy, jz, 3)*x(jx, jy, jz, 7)) &
                                /xm(jx, jy, jz, 1)**2
        end do
        end do
        end do

        do jz = 2, nz
        do jy = 2, ny
        do jx = 2, nx
            pp1 = w0(jx, jy, jz, 1)*(pr(jx + 1, jy, jz) - pr(jx - 1, jy, jz))/(2.*dx)
            pp2 = w0(jx, jy, jz, 2)*(pr(jx, jy + 1, jz) - pr(jx, jy - 1, jz))/(2.*dy)
            pp3 = w0(jx, jy, jz, 3)*(pr(jx, jy, jz + 1) - pr(jx, jy, jz - 1))/(2.*dz)
            pp4 = (pr(jx + 1, jy, jz) - pr(jx - 1, jy, jz))/(2.*dx) &
                  *(x(jx, jy, jz, 7)*(xm(jx, jy + 1, jz, 1) - xm(jx, jy - 1, jz, 1))/(2.*dy) &
                    - x(jx, jy, jz, 6)*(xm(jx, jy, jz + 1, 1) - xm(jx, jy, jz - 1, 1))/(2.*dz))
            pp5 = (pr(jx, jy + 1, jz) - pr(jx, jy - 1, jz))/(2.*dy) &
                  *(x(jx, jy, jz, 5)*(xm(jx, jy, jz + 1, 1) - xm(jx, jy, jz - 1, 1))/(2.*dz) &
                    - x(jx, jy, jz, 7)*(xm(jx + 1, jy, jz, 1) - xm(jx - 1, jy, jz, 1))/(2.*dx))
            pp6 = (pr(jx, jy, jz + 1) - pr(jx, jy, jz - 1))/(2.*dz) &
                  *(x(jx, jy, jz, 6)*(xm(jx + 1, jy, jz, 1) - xm(jx - 1, jy, jz, 1))/(2.*dx) &
                    - x(jx, jy, jz, 5)*(xm(jx, jy + 1, jz, 1) - xm(jx, jy - 1, jz, 1))/(2.*dy))
            facp(jx, jy, jz) = -(pp1 + pp2 + pp3)/xm(jx, jy, jz, 1)**2 &
                               + 2.*(pp4 + pp5 + pp6)/xm(jx, jy, jz, 1)**3
        end do
        end do
        end do

        do jz = 2, nz
        do jy = 2, ny
        do jx = 2, nx
            fact(jx, jy, jz) = x(jx, jy, jz, 5)*(xm(jx + 1, jy, jz, 2) - xm(jx - 1, jy, jz, 2))/(2.*dx) &
                               + x(jx, jy, jz, 6)*(xm(jx, jy + 1, jz, 2) - xm(jx, jy - 1, jz, 2))/(2.*dy) &
                               + x(jx, jy, jz, 7)*(xm(jx, jy, jz + 1, 2) - xm(jx, jy, jz - 1, 2))/(2.*dz)
            faci(jx, jy, jz) = fact(jx, jy, jz) - facp(jx, jy, jz)
        end do
        end do
        end do

        call bndry1(fact, 1)
        call bndry1(facp, 1)
        call bndry1(faci, 1)
        call vorticity(x)

        output = 'fac'//cn(nst)
        open (unit=8, file=output, status="unknown", form="formatted")
        write (8, 9) (((fact(jx, jy, jz), facp(jx, jy, jz), faci(jx, jy, jz), &
                        pr(jx, jy, jz), (w0(jx, jy, jz, m), m=1, 3), jx=1, mx), jy=1, my) &
                      , jz=1, mz)
9       format(7(1x, e10.4))
        close (8)
        return
    end subroutine

    subroutine energy
        implicit none
        real(kind=8) :: wyz(my, mz), fyz(my, mz), gyz(my, mz), hyz(my, mz)
        real(kind=8) :: wz(mz), fz(mz), gz(mz), hz(mz), cr1(mz), cr2(mz), crj
        real(kind=8) :: nk1(mz), nk2(mz)
        real(kind=8) :: aa, bb, cc, dd
        real(kind=8) :: wt, ft, gt, ht
        integer :: jx, jy, jz, m, k

        !  define statement functions
        !  d2fc= d2 f / dx2   with central difference
        !      d2fc(fm,f0,fp,xm1,x0,xp1)=
        !     1 2.*((fp-f0)/(xp1-x0)-(f0-fm)/(x0-xm1))/(xp1-xm1)
        !  d1fc= d f / dx  with  central difference

        !d1fc(fm, f0, fp, xm1, x0, xp1) = ((xm1 - x0)/(xp1 - x0)*(fp - f0) - (xp1 - x0)/(xm1 - x0)*(fm - f0))/(xm1 - xp1)

        do jz = 2, mz - 1
        do jy = 2, my - 1
        do jx = 2, mx - 1
            fs(jx, jy, jz) = -d1fc(pr(jx - 1, jy, jz), pr(jx, jy, jz), pr(jx + 1, jy, jz) &
                                   , xx(jx - 1), xx(jx), xx(jx + 1)) &
                             + w0(jx, jy, jz, 2)*x(jx, jy, jz, 7) - w0(jx, jy, jz, 3)*x(jx, jy, jz, 6)
            gs(jx, jy, jz) = -d1fc(pr(jx, jy - 1, jz), pr(jx, jy, jz), pr(jx, jy + 1, jz) &
                                   , yy(jy - 1), yy(jy), yy(jy + 1)) &
                             + w0(jx, jy, jz, 3)*x(jx, jy, jz, 5) - w0(jx, jy, jz, 1)*x(jx, jy, jz, 7)
            hs(jx, jy, jz) = -d1fc(pr(jx, jy, jz - 1), pr(jx, jy, jz), pr(jx, jy, jz + 1) &
                                   , zz(jz - 1), zz(jz), zz(jz + 1)) &
                             + w0(jx, jy, jz, 1)*x(jx, jy, jz, 6) - w0(jx, jy, jz, 2)*x(jx, jy, jz, 5)
        end do
        end do
        end do

        call bndry1(fs, 0)
        call bndry1(gs, 0)
        call bndry1(hs, 0)

        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            xm(jx, jy, jz, 1) = fs(jx, jy, jz)**2 + gs(jx, jy, jz)**2 + hs(jx, jy, jz)**2
            fs(jx, jy, jz) = 0.5*(x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 &
                                  + x(jx, jy, jz, 7)**2)
            gs(jx, jy, jz) = 0.5*(x(jx, jy, jz, 2)**2 + x(jx, jy, jz, 3)**2 &
                                  + x(jx, jy, jz, 4)**2)/x(jx, jy, jz, 1)
            hs(jx, jy, jz) = pr(jx, jy, jz)/(gamma - 1.)
            xm(jx, jy, jz, 2) = (w0(jx, jy, jz, 1)*x(jx, jy, jz, 5) + &
                                 w0(jx, jy, jz, 2)*x(jx, jy, jz, 6) + w0(jx, jy, jz, 3)* &
                                 x(jx, jy, jz, 7))/sqrt(2.*fs(jx, jy, jz))
        end do
        end do
        end do

        do jz = 1, mz, mz/10
            cr1(jz) = 0.
            cr2(jz) = 0.
            nk1(jz) = 1
            nk2(jz) = 1
            do jy = 1, my
            do jx = 1, mx
                crj = xm(jx, jy, jz, 2)
                if (crj .ge. 0.) then
                    nk1(jz) = nk1(jz) + 1
                    cr1(jz) = cr1(jz) + crj
                else
                    nk2(jz) = nk2(jz) + 1
                    cr2(jz) = cr2(jz) + crj
                end if
            end do
            end do
        end do

        do jz = 1, mz
        do jy = 1, my
            call integ(xm(1, jy, jz, 1), aa, xx, mx)
            call integ(fs(1, jy, jz), bb, xx, mx)
            call integ(gs(1, jy, jz), cc, xx, mx)
            call integ(hs(1, jy, jz), dd, xx, mx)
            wyz(jy, jz) = aa
            fyz(jy, jz) = bb
            gyz(jy, jz) = cc
            hyz(jy, jz) = dd
        end do
        end do

        do jz = 1, mz
            call integ(wyz(1, jz), aa, yy, my)
            call integ(fyz(1, jz), bb, yy, my)
            call integ(gyz(1, jz), cc, yy, my)
            call integ(hyz(1, jz), dd, yy, my)
            wz(jz) = aa
            fz(jz) = bb
            gz(jz) = cc
            hz(jz) = dd
        end do

        call integ(wz, wt, zz, mz)
        call integ(fz, ft, zz, mz)
        call integ(gz, gt, zz, mz)
        call integ(hz, ht, zz, mz)

        open (unit=11, file='energy.dat', status='unknown', form='formatted')
        !      write(11,9)('Force$','ME$','KE$','TE$','Time$')
        write (11, 6) wt, ft, gt, ht, time
6       format(5(1x, e13.3))

        open (unit=12, file='facur.dat', status='unknown', form='formatted')
        write (12, 7) (time)
        write (12, 11) (nk1(jz), nk2(jz), jz=1, mz, mz/10)
        write (12, 8) (cr1(jz), cr2(jz), jz=1, mz, mz/10)
7       format(f9.3)
8       format(8(1x, e9.3))
11      format(8(1x, e9.3))

        return
    end subroutine

    subroutine integ(fin, fout, x, mx)
        implicit none
        real(kind=8) :: fin(mx), fout, x(mx)
        integer :: mx

        integer :: jx

        fout = 0
        do jx = 2, mx
            fout = fout + (fin(jx - 1) + fin(jx))*(x(jx) - x(jx - 1))/2.
        end do
        return
    end subroutine

    subroutine setdt

        implicit none

        real(kind=8) ::  temp(mx,my,mz), dtmin, dxyz, test
        integer :: jx, jy, jz, m, k

        !call foreta(time, 1)
        call pressure(x, 1)
        dtmin = 1000.
        dxyz = .5*dx*dy*dz/sqrt((dx**2*dy**2 + dy**2*dz**2 + dx**2*dz**2))

        temp = dxyz/( sqrt(x(:,:,:,2)*x(:,:,:,2)+x(:,:,:,3)*x(:,:,:,3)+x(:,:,:,4)*x(:,:,:,4))/x(:,:,:,1) &
                + sqrt((x(:,:,:,5)*x(:,:,:,5)+x(:,:,:,6)*x(:,:,:,6)+x(:,:,:,7)*x(:,:,:,7)+gamma*pr)/x(:,:,:,1)) )
        
        dtmin = minval(minval(minval(temp,3),2),1)

        dt = 0.5*dtmin
        return
    end subroutine

    subroutine bndry(x, nlt)
        
        !------------------------------
        ! Set the boundaries for X's
        !------------------------------
        
        implicit none

        real(kind=8) ::  x(mx, my, mz, 8)

        integer :: jx, jy, jz, m, k
        integer :: nlt

        ! magnetosheath and -pause b.!.
        ! inflow boundary

        x(1, 2:ny, 2:nz, :) = x(2, 2:ny, 2:nz, :)
        x(mx, 2:ny, 2:nz, :) = x(nx, 2:ny, 2:nz, :)

        if (halfx) then
            do jz = 2, nz
            do jy = 2, ny
                x(mx, jy, jz, 1) = x(nx - 1, my - jy + 1, jz, 1)
                x(mx, jy, jz, 2) = -x(nx - 1, my - jy + 1, jz, 2)
                x(mx, jy, jz, 3) = -x(nx - 1, my - jy + 1, jz, 3)
                x(mx, jy, jz, 4) = x(nx - 1, my - jy + 1, jz, 4)
                x(mx, jy, jz, 5) = x(nx - 1, my - jy + 1, jz, 5)
                x(mx, jy, jz, 6) = x(nx - 1, my - jy + 1, jz, 6)
                x(mx, jy, jz, 7) = -x(nx - 1, my - jy + 1, jz, 7)
                x(mx, jy, jz, 8) = x(nx - 1, my - jy + 1, jz, 8)
            end do
            end do
        end if
        
        ! out flowing B.C.
        
        
        if (periody) then
            do m = 1, 8
            do jz = 2, nz
            do jx = 1, mx
                x(jx, 1, jz, m) = x(jx, ny, jz, m)
                x(jx, my, jz, m) = x(jx, 2, jz, m)
            end do
            end do
            end do
        else
            x(:, 1, 2:nz, :) = x(:, 2, 2:nz, :)
            x(:, my, 2:nz, :) = x(:, ny, 2:nz, :)
        end if

        x(:, :, 1, :) = x(:, :, 2, :)
        x(:, :, mz, :) = x(:, :, nz, :)

        if (halfz) then
            do jy = 1, my
            do jx = 1, mx
                x(jx, jy, mz, 1) = x(jx, jy, nz - 1, 1)
                x(jx, jy, mz, 2) = x(jx, jy, nz - 1, 2)
                x(jx, jy, mz, 3) = x(jx, jy, nz - 1, 3)
                x(jx, jy, mz, 4) = -x(jx, jy, nz - 1, 4)
                x(jx, jy, mz, 5) = -x(jx, jy, nz - 1, 5)
                x(jx, jy, mz, 6) = -x(jx, jy, nz - 1, 6)
                x(jx, jy, mz, 7) = x(jx, jy, nz - 1, 7)
                x(jx, jy, mz, 8) = x(jx, jy, nz - 1, 8)
            end do
            end do
        end if

        return
    end subroutine

    subroutine bndry1(x, nlt)
        
        !------------------------------
        ! Set the boundaries for X's
        !------------------------------

        implicit none

        real(kind=8) ::  x(mx, my, mz)
        integer :: nlt

        integer :: jx, jy, jz, m, k

        if (nlt .eq. 0) then
            ! magnetosheath and -pause b.!.

            do jz = 2, nz
            do jy = 2, ny
                x(1, jy, jz) = x(2, jy, jz)
                x(mx, jy, jz) = x(nx, jy, jz)
            end do
            end do

            if (halfx) then
                do jz = 2, nz
                do jy = 2, ny
                    x(mx, jy, jz) = x(nx - 1, my - jy + 1, jz)
                end do
                end do
            end if

            if (periody) then
                do jz = 2, nz
                do jx = 1, mx
                    x(jx, 1, jz) = x(jx, ny, jz)
                    x(jx, my, jz) = x(jx, 2, jz)
                end do
                end do
            else
                do jz = 2, nz
                do jx = 1, mx
                    x(jx, 1, jz) = x(jx, 2, jz)
                    x(jx, my, jz) = x(jx, ny, jz)
                end do
                end do
            end if

            do jy = 1, my
            do jx = 1, mx
                x(jx, jy, 1) = x(jx, jy, 2)
                if (halfz) then
                    x(jx, jy, mz) = x(jx, jy, nz - 1)
                else
                    x(jx, jy, mz) = x(jx, jy, nz)
                end if
            end do
            end do

        else
            ! magnetosheath and -pause b.!.
            do jz = 2, nz
            do jy = 2, ny
                x(1, jy, jz) = x(2, jy, jz)
                x(mx, jy, jz) = x(nx, jy, jz)
            end do
            end do

            if (halfx) then
                do jz = 2, nz
                do jy = 2, ny
                    x(mx, jy, jz) = -x(nx - 1, my - jy + 1, jz)
                end do
                end do
            end if

            if (periody) then
                do jz = 2, nz
                do jx = 1, mx
                    x(jx, 1, jz) = x(jx, ny, jz)
                    x(jx, my, jz) = x(jx, 2, jz)
                end do
                end do
            else
                do jz = 2, nz
                do jx = 1, mx
                    x(jx, 1, jz) = x(jx, 2, jz)
                    x(jx, my, jz) = x(jx, ny, jz)
                end do
                end do
            end if

            do jy = 1, my
            do jx = 1, mx
                x(jx, jy, 1) = x(jx, jy, 2)
                if (halfz) then
                    x(jx, jy, mz) = -x(jx, jy, nz - 1)
                else
                    x(jx, jy, mz) = x(jx, jy, nz)
                end if
            end do
            end do
        end if

        return
    end subroutine

    subroutine flux(x, fs, gs, hs, nnx, nny, nnz, m, mm)
        
        !-------------------------------
        !  Calculate fluxes
        !  Notations: X1    X2     X3     X4     X5  X6  x7
        !             rho   rhovx  rhovy  rhovz  bx  by  bz  e
        !-------------------------------
        

        implicit none

        real(kind=8) :: x(mx, my, mz, 8), fs(mx, my, mz)
        real(kind=8) :: hs(mx, my, mz), gs(mx, my, mz)
        integer :: nnx, nny, nnz, m, mm

        real(kind=8) :: vcrbz, vcrby, vcrbx
        real(kind=8) :: b2, bdotv, eng
        integer :: jx, jy, jz

        if (m .eq. 1) then
            ! [1] Continuity eq.
            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                fs(jx, jy, jz) = x(jx, jy, jz, 2)
                gs(jx, jy, jz) = x(jx, jy, jz, 3)
                hs(jx, jy, jz) = x(jx, jy, jz, 4)
            end do
            end do
            end do
            !$OMP END PARALLEL DO
            
        else if (m .eq. 2) then
            ![2] Momentum eq.
            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                fs(jx, jy, jz) = x(jx, jy, jz, 2)**2/x(jx, jy, jz, 1) + pr(jx, jy, jz) &
                                 + .5*b2 - x(jx, jy, jz, 5)**2
                gs(jx, jy, jz) = x(jx, jy, jz, 2)*x(jx, jy, jz, 3)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 6)*x(jx, jy, jz, 5)
                hs(jx, jy, jz) = x(jx, jy, jz, 2)*x(jx, jy, jz, 4)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 7)*x(jx, jy, jz, 5)
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 3) then
            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                fs(jx, jy, jz) = x(jx, jy, jz, 2)*x(jx, jy, jz, 3)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 6)*x(jx, jy, jz, 5)
                gs(jx, jy, jz) = x(jx, jy, jz, 3)**2/x(jx, jy, jz, 1) + pr(jx, jy, jz) &
                                 + .5*b2 - x(jx, jy, jz, 6)**2
                hs(jx, jy, jz) = x(jx, jy, jz, 4)*x(jx, jy, jz, 3)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 6)*x(jx, jy, jz, 7)
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 4) then
            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx

                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                fs(jx, jy, jz) = x(jx, jy, jz, 2)*x(jx, jy, jz, 4)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 5)*x(jx, jy, jz, 7)
                gs(jx, jy, jz) = x(jx, jy, jz, 3)*x(jx, jy, jz, 4)/x(jx, jy, jz, 1) &
                                 - x(jx, jy, jz, 6)*x(jx, jy, jz, 7)
                hs(jx, jy, jz) = x(jx, jy, jz, 4)**2/x(jx, jy, jz, 1) + pr(jx, jy, jz) &
                                 + .5*b2 - x(jx, jy, jz, 7)**2
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 5) then
            ! [3] Magnetic induction eq.

            !call current(x, mm)
            !call foreta(time, mm)

            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                vcrbz = ((x(jx, jy, jz, 2) - di*w0(jx, jy, jz, 1))*x(jx, jy, jz, 6) &
                         - (x(jx, jy, jz, 3) - di*w0(jx, jy, jz, 2))*x(jx, jy, jz, 5)) &
                        /x(jx, jy, jz, 1)
                vcrby = ((x(jx, jy, jz, 4) - di*w0(jx, jy, jz, 3))*x(jx, jy, jz, 5) &
                         - (x(jx, jy, jz, 2) - di*w0(jx, jy, jz, 1))*x(jx, jy, jz, 7)) &
                        /x(jx, jy, jz, 1)

                fs(jx, jy, jz) = 0.
                gs(jx, jy, jz) = -vcrbz + etaf(jx, jy, jz)*w0(jx, jy, jz, 3)
                hs(jx, jy, jz) = vcrby - etaf(jx, jy, jz)*w0(jx, jy, jz, 2)
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 6) then
            ! [3] Magnetic induction eq.

            !call current(x, mm)
            !call foreta(time, mm)

            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx

                

                vcrbz = ((x(jx, jy, jz, 2) - di*w0(jx, jy, jz, 1))*x(jx, jy, jz, 6) &
                         - (x(jx, jy, jz, 3) - di*w0(jx, jy, jz, 2))*x(jx, jy, jz, 5)) &
                        /x(jx, jy, jz, 1)
                vcrbx = ((x(jx, jy, jz, 3) - di*w0(jx, jy, jz, 2))*x(jx, jy, jz, 7) &
                         - (x(jx, jy, jz, 4) - di*w0(jx, jy, jz, 3))*x(jx, jy, jz, 6)) &
                        /x(jx, jy, jz, 1)

                fs(jx, jy, jz) = vcrbz - etaf(jx, jy, jz)*w0(jx, jy, jz, 3)
                gs(jx, jy, jz) = 0.
                hs(jx, jy, jz) = -vcrbx + etaf(jx, jy, jz)*w0(jx, jy, jz, 1)
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 7) then
            ! [3] Magnetic induction eq.

            !call current(x, mm)
            !call foreta(time, mm)

            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                vcrby = ((x(jx, jy, jz, 4) - di*w0(jx, jy, jz, 3))*x(jx, jy, jz, 5) &
                         - (x(jx, jy, jz, 2) - di*w0(jx, jy, jz, 1))*x(jx, jy, jz, 7)) &
                        /x(jx, jy, jz, 1)
                vcrbx = ((x(jx, jy, jz, 3) - di*w0(jx, jy, jz, 2))*x(jx, jy, jz, 7) &
                         - (x(jx, jy, jz, 4) - di*w0(jx, jy, jz, 3))*x(jx, jy, jz, 6)) &
                        /x(jx, jy, jz, 1)
                fs(jx, jy, jz) = -vcrby + etaf(jx, jy, jz)*w0(jx, jy, jz, 2)
                gs(jx, jy, jz) = vcrbx - etaf(jx, jy, jz)*w0(jx, jy, jz, 1)
                hs(jx, jy, jz) = 0.
            end do
            end do
            end do
            !$OMP END PARALLEL DO

        else if (m .eq. 8) then
            ! [4] Energy eq.

            !$OMP PARALLEL DO
            do jz = 1, nnz
            do jy = 1, nny
            do jx = 1, nnx
                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                bdotv = (x(jx, jy, jz, 5)*x(jx, jy, jz, 2) + x(jx, jy, jz, 6)*x(jx, jy, jz, 3) &
                         + x(jx, jy, jz, 7)*x(jx, jy, jz, 4))/x(jx, jy, jz, 1)
                eng = x(jx, jy, jz, 8) + pr(jx, jy, jz) + .5*b2

                fs(jx, jy, jz) = eng*x(jx, jy, jz, 2)/x(jx, jy, jz, 1) &
                                 - bdotv*x(jx, jy, jz, 5) &
                                 + etaf(jx, jy, jz)*(w0(jx, jy, jz, 2)*x(jx, jy, jz, 7) &
                                                     - w0(jx, jy, jz, 3)*x(jx, jy, jz, 6))
                gs(jx, jy, jz) = eng*x(jx, jy, jz, 3)/x(jx, jy, jz, 1) &
                                 - bdotv*x(jx, jy, jz, 6) &
                                 + etaf(jx, jy, jz)*(w0(jx, jy, jz, 3)*x(jx, jy, jz, 5) &
                                                     - w0(jx, jy, jz, 1)*x(jx, jy, jz, 7))
                hs(jx, jy, jz) = eng*x(jx, jy, jz, 4)/x(jx, jy, jz, 1) &
                                 - bdotv*x(jx, jy, jz, 7) &
                                 + etaf(jx, jy, jz)*(w0(jx, jy, jz, 1)*x(jx, jy, jz, 6) &
                                                     - w0(jx, jy, jz, 2)*x(jx, jy, jz, 5))

            end do
            end do
            end do
            !$OMP END PARALLEL DO

        end if

        return
    end subroutine

    subroutine current(x, mm)
        implicit none

        real(kind=8) :: x(mx, my, mz, 8)
        integer :: mm

        integer :: m, jx, jy, jz

        if (mm .eq. 1) then
            do jz = 2, nz
            do jy = 2, ny
            do jx = 2, nx
                w0(jx, jy, jz, 1) = .5*(x(jx, jy + 1, jz, 7) - x(jx, jy - 1, jz, 7))/dy &
                                    - .5*(x(jx, jy, jz + 1, 6) - x(jx, jy, jz - 1, 6))/dz
                w0(jx, jy, jz, 2) = -.5*(x(jx + 1, jy, jz, 7) - x(jx - 1, jy, jz, 7))/dx &
                                    + .5*(x(jx, jy, jz + 1, 5) - x(jx, jy, jz - 1, 5))/dz
                w0(jx, jy, jz, 3) = .5*(x(jx + 1, jy, jz, 6) - x(jx - 1, jy, jz, 6))/dx &
                                    - .5*(x(jx, jy + 1, jz, 5) - x(jx, jy - 1, jz, 5))/dy
            end do
            end do
            end do
            
            ! boundary at jx=1,mx

            do m = 1, 3
            do jz = 2, nz
            do jy = 2, ny
                w0(1, jy, jz, m) = 2.*w0(2, jy, jz, m) - w0(3, jy, jz, m)
                w0(mx, jy, jz, m) = 2.*w0(nx, jy, jz, m) - w0(nx - 1, jy, jz, m)
            end do
            end do
            end do

            if (halfx) then
                do jz = 2, nz
                do jy = 2, ny
                    w0(mx, jy, jz, 1) = w0(nx - 1, my - jy + 1, jz, 1)
                    w0(mx, jy, jz, 2) = w0(nx - 1, my - jy + 1, jz, 2)
                    w0(mx, jy, jz, 3) = -w0(nx - 1, my - jy + 1, jz, 3)
                end do
                end do
            end if

            ! boundary at jy=1,my
            if (periody) then
                do m = 1, 3
                do jz = 2, nz
                do jx = 1, mx
                    w0(jx, 1, jz, m) = w0(jx, ny, jz, m)
                    w0(jx, my, jz, m) = w0(jx, 2, jz, m)
                end do
                end do
                end do
            else
                do m = 1, 3
                do jz = 2, nz
                do jx = 1, mx
                    w0(jx, 1, jz, m) = w0(jx, 2, jz, m)
                    w0(jx, my, jz, m) = w0(jx, ny, jz, m)
                end do
                end do
                end do
            end if

            ! b.!. at jz=1
            do m = 1, 3
            do jy = 1, my
            do jx = 1, mx
                w0(jx, jy, 1, m) = w0(jx, jy, 2, m)
                w0(jx, jy, mz, m) = w0(jx, jy, nz, m)
            end do
            end do
            end do

            ! b.!. at jz=mz
            if (halfz) then
                do jy = 1, my
                do jx = 1, mx
                    w0(jx, jy, mz, 1) = w0(jx, jy, nz - 1, 1)
                    w0(jx, jy, mz, 2) = w0(jx, jy, nz - 1, 2)
                    w0(jx, jy, mz, 3) = -w0(jx, jy, nz - 1, 3)
                end do
                end do
            end if

        else
            do jz = 2, nz - 1
            do jy = 2, ny - 1
            do jx = 2, nx - 1
                w0(jx, jy, jz, 1) = .5*(x(jx, jy + 1, jz, 7) - x(jx, jy - 1, jz, 7))/dy &
                                    - .5*(x(jx, jy, jz + 1, 6) - x(jx, jy, jz - 1, 6))/dz
                w0(jx, jy, jz, 2) = -.5*(x(jx + 1, jy, jz, 7) - x(jx - 1, jy, jz, 7))/dx &
                                    + .5*(x(jx, jy, jz + 1, 5) - x(jx, jy, jz - 1, 5))/dz
                w0(jx, jy, jz, 3) = .5*(x(jx + 1, jy, jz, 6) - x(jx - 1, jy, jz, 6))/dx &
                                    - .5*(x(jx, jy + 1, jz, 5) - x(jx, jy - 1, jz, 5))/dy
            end do
            end do
            end do

            ! boundary at jx=1,nx
            do m = 1, 3
            do jz = 2, nz - 1
            do jy = 2, ny - 1
                w0(1, jy, jz, m) = w0(2, jy, jz, m)
                w0(nx, jy, jz, m) = w0(nx - 1, jy, jz, m)
            end do
            end do
            end do

            if (halfx) then
                do jz = 2, nz
                do jy = 2, ny
                    w0(nx, jy, jz, 1) = w0(nx - 1, ny - jy + 1, jz, 1)
                    w0(nx, jy, jz, 2) = w0(nx - 1, ny - jy + 1, jz, 2)
                    w0(nx, jy, jz, 3) = -w0(nx - 1, ny - jy + 1, jz, 3)
                end do
                end do
            end if

            ! boundary at jy=1,ny
            if (periody) then
                do m = 1, 3
                do jz = 2, nz - 1
                do jx = 1, nx
                    w0(jx, 1, jz, m) = w0(jx, ny - 1, jz, m)
                    w0(jx, ny, jz, m) = w0(jx, 2, jz, m)
                end do
                end do
                end do
            else
                do m = 1, 3
                do jz = 2, nz - 1
                do jx = 1, nx
                    w0(jx, 1, jz, m) = w0(jx, 2, jz, m)
                    w0(jx, ny, jz, m) = w0(jx, ny - 1, jz, m)
                end do
                end do
                end do
            end if

            do m = 1, 3
            do jy = 1, ny
            do jx = 1, nx
                w0(jx, jy, 1, m) = w0(jx, jy, 2, m)
                w0(jx, jy, nz, m) = w0(jx, jy, nz - 1, m)
            end do
            end do
            end do

            ! b.!. at jz=nz
            if (halfz) then
                do jy = 1, ny
                do jx = 1, nx
                    w0(jx, jy, nz, 1) = w0(jx, jy, nz - 1, 1)
                    w0(jx, jy, nz, 2) = w0(jx, jy, nz - 1, 2)
                    w0(jx, jy, nz, 3) = -w0(jx, jy, nz - 1, 3)
                end do
                end do
            end if

        end if

        return
    end subroutine

    subroutine vorticity(x)

        implicit none

        real(kind=8) ::  x(mx, my, mz, 8)
        
        real(kind=8) :: velocity(mx, my, mz, 3), vor(mx, my, mz, 3)

        integer :: jx, jy, jz, m

        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            velocity(jx, jy, jz, 1) = x(jx, jy, jz, 2)/x(jx, jy, jz, 1)
            velocity(jx, jy, jz, 2) = x(jx, jy, jz, 3)/x(jx, jy, jz, 1)
            velocity(jx, jy, jz, 3) = x(jx, jy, jz, 4)/x(jx, jy, jz, 1)
        end do
        end do
        end do

        do jz = 2, nz
        do jy = 2, ny
        do jx = 2, nx
            vor(jx, jy, jz, 1) = .5*(velocity(jx, jy + 1, jz, 3) - velocity(jx, jy - 1, jz, 3))/dy &
                                - .5*(velocity(jx, jy, jz + 1, 2) - velocity(jx, jy, jz - 1, 2))/dz
            vor(jx, jy, jz, 2) = -.5*(velocity(jx + 1, jy, jz, 3) - velocity(jx - 1, jy, jz, 3))/dx &
                                + .5*(velocity(jx, jy, jz + 1, 1) - velocity(jx, jy, jz - 1, 1))/dz
            vor(jx, jy, jz, 3) = .5*(velocity(jx + 1, jy, jz, 2) - velocity(jx - 1, jy, jz, 2))/dx &
                                - .5*(velocity(jx, jy + 1, jz, 1) - velocity(jx, jy - 1, jz, 1))/dy
        end do
        end do
        end do

        ! boundary at jx=1,mx
        do m = 1, 3
        do jz = 2, nz
        do jy = 2, ny
            vor(1, jy, jz, m) = 2.*vor(2, jy, jz, m) - vor(3, jy, jz, m)
            vor(mx, jy, jz, m) = 2.*vor(nx, jy, jz, m) - vor(nx - 1, jy, jz, m)
        end do
        end do
        end do

        if (halfx) then
            do jz = 2, nz
            do jy = 2, ny
                vor(mx, jy, jz, 1) = vor(nx - 1, my - jy + 1, jz, 1)
                vor(mx, jy, jz, 2) = vor(nx - 1, my - jy + 1, jz, 2)
                vor(mx, jy, jz, 3) = -vor(nx - 1, my - jy + 1, jz, 3)
            end do
            end do
        end if

        ! boundary at jy=1,my
        if (periody) then
            do m = 1, 3
            do jz = 2, nz
            do jx = 1, mx
                vor(jx, 1, jz, m) = vor(jx, ny, jz, m)
                vor(jx, my, jz, m) = vor(jx, 2, jz, m)
            end do
            end do
            end do
        else
            do m = 1, 3
            do jz = 2, nz
            do jx = 1, mx
                vor(jx, 1, jz, m) = vor(jx, 2, jz, m)
                vor(jx, my, jz, m) = vor(jx, ny, jz, m)
            end do
            end do
            end do
        end if

        ! b.!. at jz=1
        do m = 1, 3
        do jy = 1, my
        do jx = 1, mx
            vor(jx, jy, 1, m) = vor(jx, jy, 2, m)
            vor(jx, jy, mz, m) = vor(jx, jy, nz, m)
        end do
        end do
        end do

        ! b.!. at jz=mz
        if (halfz) then
            do jy = 1, my
            do jx = 1, mx
                vor(jx, jy, mz, 1) = -vor(jx, jy, nz - 1, 1)
                vor(jx, jy, mz, 2) = -vor(jx, jy, nz - 1, 2)
                vor(jx, jy, mz, 3) = vor(jx, jy, nz - 1, 3)
            end do
            end do
        end if

        return
    end subroutine

    !      subroutine readin
    !      include 'ma3ds1.for'
    !      include 'ma3ds2.for'
    !      character*8 contin
    !    character*3 cn

    !      contin='m3d'//cn(nst)
    !     open(unit=8,file=contin,status="unknown",form="unformatted")
    !    open(unit=8,file='continue',status="unknown",form="unformatted")
    !      read(8)ncase,nstep,time,nst
    !      read(8)x
    !     close(8)
    !    nst=1
    !      return
    !      end
    !
    !cyg-------------------------------
    subroutine readin(nst3, cont2, dtime1)

        implicit none

        integer :: nst3, cont2
        real(kind=8) :: dtime1

        character*8 contin
        integer :: nst

        nst = nst3

        if (cont2 .eq. 1) then
            open (unit=8, file='continue', status="unknown", form="unformatted")
            read (8) ncase, nstep, time, nst
            read (8) x
            close (8)
        else
            contin = 'm3d'//cn(nst)
            open (unit=8, file=contin, status="unknown", form="unformatted")
            read (8) ncase, nstep, time
            read (8) x
            nst = ceiling(time)/dtime1 + 1
            close (8)
        end if

        return
    end subroutine

    !cyg----------------------
    subroutine recrd

        implicit none

        character*8 output

        output = 'm3d'//cn(nst)
        open (unit=8, file=output, status="unknown", form="unformatted")
        write (8) ncase, nstep, time
        write (8) x
        close (8)
        return
    end subroutine

    subroutine recrd1

        implicit none

        character*8 output
        integer :: jx, jy, jz, m

        output = 'm3ds'//cn(nst)
        !call current(x, 1)
        open (unit=8, file=output, status="unknown", form="formatted")
        write (8, 9) ((((x(jx, jy, jz, m), m=1, 8), (w0(jx, jy, jz, m), m=1, 3), jx=1, mx), jy=1, my), jz=1, mz)
9       format(11(1x, e10.4))
        close (8)
        return
    end subroutine

    subroutine gridpnt

        implicit none

        integer :: jx, jy, jz, m

        nnxh = (mx + 1)/2
        nnx = 2*(nnxh - 1)
        nnxq = nnxh/2
        nnyh = (my + 1)/2
        nny = 2*(nnyh - 1)
        nnyq = nnyh/2
        nnzh = (mz + 1)/2
        nnz = 2*(nnzh - 1)
        nnzq = nnzh/2

        dx = (xmax - xmin)/float(nnx)

        do jx = 1, mx
            xx(jx) = xmin + (jx - 1)*dx
        end do

        dy = (ymax - ymin)/float(nny)

        do jy = 1, my
            yy(jy) = ymin + (jy - 1)*dy
        end do

        dz = (zmax - zmin)/float(nnz)

        do jz = 1, mz
            zz(jz) = zmin + (jz - 1)*dz
        end do

        open (unit=11, file='grid.dat', status='unknown', form='formatted')
        write (11, 99) (xx(jx), jx=1, mx), (zz(jz), jz=1, mz)
99      format(5(1x, e10.4))

        return
    end subroutine

    subroutine smthxyz(x, weight, num)

        implicit none
        real(kind=8) ::  x(mx, my, mz, 8)
        real(kind=8) :: weight   ! weight of the central element
        integer :: num  ! number of the smooth operation

        integer :: jx, jy, jz, m, k
        real(kind=8) :: average

        do k = 1, num
            do m = 1, 8
                do jz = 2, nz
                do jy = 2, ny
                do jx = 2, nx
                    ! differ(jx, jy, jz) = ((x(jx + 1, jy, jz, m) + x(jx - 1, jy, jz, m)) &
                    !                      + (x(jx, jy + 1, jz, m) + x(jx, jy - 1, jz, m)) &
                    !                      + (x(jx, jy, jz + 1, m) + x(jx, jy, jz - 1, m)) - 6.*x(jx, jy, jz, m))
                    ! x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./48.0)*differ(jx, jy, jz)

                    average = ((x(jx + 1, jy, jz, m) + x(jx - 1, jy, jz, m)) &
                                + (x(jx, jy + 1, jz, m) + x(jx, jy - 1, jz, m)) &
                                + (x(jx, jy, jz + 1, m) + x(jx, jy, jz - 1, m)))/6.0
                    x(jx, jy, jz, m) = weight*x(jx, jy, jz, m) + (1 - weight)*average
                end do
                end do
                end do

                ! do jz = 2, nz
                ! do jy = 2, ny
                ! do jx = 2, nsmthx
                !     theta = 3.1415926*(jx - 2)/(nsmthx - 3)
                !     x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./96.0)*(.5*(1.+cos(theta)))*differ(jx, jy, jz, 1)
                ! end do
                ! end do
                ! end do

                ! if (.not. halfx) then
                !     do jz = 2, nz
                !     do jy = 2, ny
                !     do jx = mx - nsmthx + 1, nx
                !         theta = 3.1415926*(mx - jx - 1)/(nsmthx - 3)
                !         x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./96.0)*(.5*(1.+cos(theta)))*differ(jx, jy, jz, 1)
                !     end do
                !     end do
                !     end do
                ! end if

                ! if (.not. periody) then
                !     do jz = 2, nz
                !     do jy = 2, nsmthy
                !     do jx = 2, nx
                !         theta = 3.1415926*(jy - 2)/(nsmthy - 3)
                !         x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./96.0)*(.5*(1.+cos(theta)))*differ(jx, jy, jz, 1)
                !     end do
                !     end do
                !     end do

                !     do jz = 2, nz
                !     do jy = my - nsmthy + 1, ny
                !     do jx = 2, nx
                !         theta = 3.1415926*(my - jy - 1)/(nsmthy - 3)
                !         x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./96.0)*(.5*(1.+cos(theta)))*differ(jx, jy, jz, 1)
                !     end do
                !     end do
                !     end do
                ! end if

                ! do jz = 2, nz
                ! do jy = 2, ny
                ! do jx = 2, nx
                !     ! theta=2.*3.1415926*zz(jz)/zmin
                !     x(jx, jy, jz, m) = x(jx, jy, jz, m) + (1./48.0)*differ(jx, jy, jz)
                !     !     1     +(1./48.0)*(2.+cos(theta))/3.*differ(jx,jy,jz,1)
                ! end do
                ! end do
                ! end do

                !      if(.not.halfz) then
                !      do 53 jz=mz-nsmthz+1,nz
                !      do 53 jy=2,ny
                !      do 53 jx=2,nx
                !      theta=3.1415926*(mz-jz-1)/(nsmthz-3)
                !      x(jx,jy,jz,m)=x(jx,jy,jz,m)
                !     1     +(1./96.0)*(.5*(1.+cos(theta)))*differ(jx,jy,jz,1)
                !   53 continue
                !      endif
            end do
            call bndry(x, 2)
        end do

        return
    end subroutine

    subroutine avrg(x, caf1)
        ! combine the avrg1 and avrg2 to avrg

        implicit none

        real(kind=8) ::  x(mx, my, mz, 8), caf1

        integer :: jx, jy, jz, m

        do m = 1, 8
            do jz = 2, nz
            do jy = 2, ny
            do jx = 2, nx
                x(jx, jy, jz, m) = caf1*x(jx, jy, jz, m) + (1 - caf1)* &
                                    ((x(jx + 1, jy, jz, m) + x(jx - 1, jy, jz, m)) &
                                     + (x(jx, jy + 1, jz, m) + x(jx, jy - 1, jz, m)) &
                                     + (x(jx, jy, jz + 1, m) + x(jx, jy, jz - 1, m)))/6.0
            end do
            end do
            end do

            ! combine the two loop un and down to one
            ! do jz = 2, nz
            ! do jy = 2, ny
            ! do jx = 2, nx
            !     x(jx, jy, jz, m) = w0(jx, jy, jz, 1)
            ! end do
            ! end do
            ! end do
        end do

        call bndry(x, 2)

        return
    end subroutine

    subroutine avrg2(x, caf1)
        ! average x( , , ,5:7)

        implicit none

        real(kind=8) :: x(mx, my, mz, 8), caf1

        integer :: jx, jy, jz, m

        do m = 5, 7
            do jz = 2, nz
            do jy = 2, ny
            do jx = 2, nx
                x(jx, jy, jz, m) = caf1*x(jx, jy, jz, m) + (1 - caf1)* &
                                    ((x(jx + 1, jy, jz, m) + x(jx - 1, jy, jz, m)) &
                                     + (x(jx, jy + 1, jz, m) + x(jx, jy - 1, jz, m)) &
                                     + (x(jx, jy, jz + 1, m) + x(jx, jy, jz - 1, m)))/6.0
            end do
            end do
            end do

            ! combine the two loop un and down to one
            ! do jz = 2, nz
            ! do jy = 2, ny
            ! do jx = 2, nx
            !     x(jx, jy, jz, m) = w0(jx, jy, jz, 1)
            ! end do
            ! end do
            ! end do
        end do

        call bndry(x, 2)

        return
    end subroutine

    subroutine smthf(x, caf1)

        implicit none

        real(kind=8) ::  x(mx, my, mz, 8), caf1

        integer :: jx, jy, jz, m
        real(kind=8) :: wh(mx, my, mz, 3)

        do m = 1, 8
            do jz = 2, nz
            do jx = 2, nx
            do jy = 2, ny
                wh(jx, jy, jz, 1) = (x(jx + 1, jy, jz, m) - x(jx, jy, jz, m))* &
                                    (x(jx, jy, jz, m) - x(jx - 1, jy, jz, m))
                wh(jx, jy, jz, 2) = (x(jx, jy + 1, jz, m) - x(jx, jy, jz, m))* &
                                    (x(jx, jy, jz, m) - x(jx, jy - 1, jz, m))
                wh(jx, jy, jz, 3) = (x(jx, jy, jz + 1, m) - x(jx, jy, jz, m))* &
                                    (x(jx, jy, jz, m) - x(jx, jy, jz - 1, m))
            end do
            end do
            end do

            call bndry1(wh(1, 1, 1, 1), 0)
            call bndry1(wh(1, 1, 1, 2), 0)
            call bndry1(wh(1, 1, 1, 3), 0)

            do jz = 2, nz
            do jx = 2, nx
            do jy = 2, ny
                if ((wh(jx, jy, jz, 1) .lt. 0 .and. &
                     (wh(jx + 1, jy, jz, 1) .lt. 0 .or. wh(jx - 1, jy, jz, 1) .lt. 0)) &
                    .or. (wh(jx, jy, jz, 2) .lt. 0 .and. &
                          (wh(jx, jy + 1, jz, 2) .lt. 0 .or. wh(jx, jy - 1, jz, 2) .lt. 0)) &
                    .or. (wh(jx, jy, jz, 3) .lt. 0 .and. &
                          (wh(jx, jy, jz + 1, 3) .lt. 0 .or. wh(jx, jy, jz - 1, 3) .lt. 0))) then
                    w0(jx, jy, jz, 1) = caf1*x(jx, jy, jz, m) + (1 - caf1)* &
                                        ((x(jx + 1, jy, jz, m) + x(jx - 1, jy, jz, m)) &
                                         + (x(jx, jy + 1, jz, m) + x(jx, jy - 1, jz, m)) &
                                         + (x(jx, jy, jz + 1, m) + x(jx, jy, jz - 1, m)))/6.0
                else
                    w0(jx, jy, jz, 1) = x(jx, jy, jz, m)
                end if
            end do
            end do
            end do

            do jz = 2, mz - 1
            do jx = 2, nx
            do jy = 2, ny
                x(jx, jy, jz, m) = w0(jx, jy, jz, 1)
            end do
            end do
            end do
        end do

        call bndry(x, 2)
        return
    end subroutine

    subroutine pressure(x, mm)

        implicit none

        real(kind=8) :: x(mx, my, mz, 8)
        integer :: mm

        real(kind=8) :: rv2, b2
        integer :: jx, jy, jz, m

        if (mm .eq. 1) then
            do jz = 1, mz
            do jy = 1, my
            do jx = 1, mx
                rv2 = (x(jx, jy, jz, 2)**2 + x(jx, jy, jz, 3)**2 + x(jx, jy, jz, 4)**2) &
                      /x(jx, jy, jz, 1)
                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                pr(jx, jy, jz) = (gamma - 1.)*(x(jx, jy, jz, 8) - .5*rv2 - .5*b2)
            end do
            end do
            end do

        else
            do jz = 1, nz
            do jy = 1, ny
            do jx = 1, nx
                rv2 = (x(jx, jy, jz, 2)**2 + x(jx, jy, jz, 3)**2 + x(jx, jy, jz, 4)**2) &
                      /x(jx, jy, jz, 1)
                b2 = x(jx, jy, jz, 5)**2 + x(jx, jy, jz, 6)**2 + x(jx, jy, jz, 7)**2
                pr(jx, jy, jz) = (gamma - 1.)*(x(jx, jy, jz, 8) - .5*rv2 - .5*b2)
            end do
            end do
            end do

        end if

        call positive(pr, 1.d-5)

        return
    end subroutine

    subroutine positive(fn, c)

        implicit none

        real(kind=8) :: fn(mx, my, mz), c

        character*15 out
        integer :: jx, jy, jz

        do jz = 1, mz
        do jy = 1, my
        do jx = 1, mx
            if (fn(jx, jy, jz) .lt. 0.0) then
                out = 'finaltime.txt'
                open (unit=8, file=out, status="unknown", form="formatted")
                write (8, 20) time
20              format(9(1x, e10.4))
                !    write(*,*)'finaltime=',time
                !    stop
            end if

            if (fn(jx, jy, jz) .lt. c) then
                fn(jx, jy, jz) = c
            end if
        end do
        end do
        end do

        return
    end subroutine

    subroutine foreta(t, mm)
        !t seems not been used

        implicit none

        real(kind=8) :: t
        integer :: mm

        real(kind=8) :: etam(my), etac, etal, etab, xlen, xtrig, ztrig, awx, awz, etax, etaz, alpha0
        integer :: jx, jy, jz

        if(.not. is_foreta) return

        etac = 0.5

        !cjx  attentiion
        etal = 0.000
        etab = 0.00005
        alpha0 = 2.0
        !cjx  attentiion
        xlen = 15.0
        xtrig = 0.0
        !cjx  attentiion
        ztrig = 0.0
        !cjx  attentiion
        awx = aw
        awz = 2.*awx

        if (mm .eq. 1) then
            do jy = 1, my
                if (abs(yy(jy)) .le. xlen) then
                    etam(jy) = 1.
                else
                    etam(jy) = (1.-tanh(abs(yy(jy)) - xlen)**2)
                end if
            end do

            do jz = 1, mz
            do jy = 1, my
            do jx = 1, mx
                ! nonlinear resistivity
                !
                ! localized resistivity perturbation
                if (abs(xx(jx) - xtrig) .le. 0.5) then
                    etax = 1.0
                else
                    etax = 1.-tanh((xx(jx) - xtrig)/awx)**2
                end if

                if (abs(abs(zz(jz)) - ztrig) .le. 1.) then
                    etaz = 1.0
                else
                    etaz = 1.-tanh((abs(zz(jz)) - ztrig)/awz)**2
                end if

                etaf(jx, jy, jz) = etab + etam(jy)*etal*etax*etaz
            end do
            end do
            end do
        else
            do jy = 1, ny
                if (abs(yy(jy) + dy/2.) .le. xlen) then
                    etam(jy) = 1.
                else
                    etam(jy) = (1.-tanh(abs(yy(jy) + dy/2.) - xlen)**2)
                end if
            end do

            do jz = 1, nz
            do jy = 1, ny
            do jx = 1, nx
                ! localized resistivity perturbation
                if (abs(xx(jx) + dx/2.-xtrig) .le. 0.5) then
                    etax = 1.0
                else
                    etax = 1.-tanh((xx(jx) + dx/2.-xtrig)/awx)**2
                end if

                if (abs(abs(zz(jz) + dz/2.) - ztrig) .le. 1.) then
                    etaz = 1.0
                else
                    etaz = 1.-tanh((zz(jz) + dz/2.-ztrig)/awz)**2
                end if

                etaf(jx, jy, jz) = etab + etam(jy)*etal*etax*etaz
            end do
            end do
            end do
        end if

        return
    end subroutine

    !cyg---------------------------------------
    subroutine incident_plasma(x, xi, t, t1, Io)
        implicit none

        real(kind=8) :: x(mx, my, mz, 8), xi(mx, my, mz, 8), t, t1
        integer :: Io

        real(kind=8) :: vs0, lsx, lsy, lsz, tao, xx0, xx1, yy0, yy1, zz0, zz1
        real(kind=8) :: theta_in, phi_in, vx0_in, vy0_in, vz0_in, vx1_in, vy1_in, vz1_in
        integer :: a, b
        integer :: jx, jy, jz, m

        vs0 = 0.6
        lsz = 5.
        lsx = 1.
        lsy = 1.
        tao = 10
        t0 = t1
        xx0 = 0.
        xx1 = xx0

        yy0 = 0.
        yy1 = -yy0

        zz0 = 0.
        zz1 = -zz0

        theta_in = 10
        theta_in = theta_in*pi/180

        phi_in = 0
        phi_in = phi_in*pi/180

        vx0_in = vs0*sin(theta_in)*cos(phi_in)
        vy0_in = vs0*sin(theta_in)*sin(phi_in)
        vz0_in = vs0*cos(theta_in)

        vx1_in = vx0_in
        vy1_in = vy0_in
        vz1_in = -vz0_in

        if (Io .eq. 0) then
            a = 1
            b = 0
        end if

        if (Io .eq. 1) then
            a = 1
            b = 0
        end if

        

            if (t < tao) then
                do jz = 1, mz
                do jx = 1, mx
                    x(jx, :, jz, 3) = x(jx, :, jz, 1)*vs0*exp(-(((xx(jx) - xx0)/lsx)**2  + ((zz(jz) - zz0)/lsz)**2)) &
                                    *tanh(t/tao)
                end do
                end do
            end if

        

        return
    end subroutine
    !cyg-------------------------

end module Custom_subroutines
