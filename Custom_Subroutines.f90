module custom_subroutines
    use define_variables
    use custom_functions
    implicit none

contains

    subroutine initialize
        implicit none

        integer :: i
        real(kind=8) :: Bh = -1.0d0 ! Harris sheet magnetic field

        ! constants and normalizing parameters
        x0 = 3.4d5
        t0 = 1.42344d-4
        n0 = 5.0d6
        v0 = x0 / t0
        E0 = m0 * v0 / (q0 * t0)
        B0 = m0 / (q0 * t0)
        pr0 = m0 * n0 * v0**2
        Tem0 = m0 * v0**2
        j0 = B0 / (mu0 * x0)
        eta0 = m0 / (q0**2 * n0 * t0)
        R0 = eta0 * q0 * n0 * j0
        const1 = mu0 * x0**2 * n0 * q0**2 / m0   ! used in calculating current density j
        const2 = t0**2 / (mu0 * eps0 * x0**2)    ! used in calculating electric field E
        const3 = R0 / (m0 * n0 * v0 / t0)   ! used in fraction term Re and Ri
        qe = -1.0d0
        qi = 1.0d0
        me = 1.0d0
        mi = 1864.0d0

        print *, "x0 = ", x0
        print *, "t0 = ", t0
        print *, "n0 = ", n0
        print *, "v0 = ", v0
        print *, "E0 = ", E0
        print *, "B0 = ", B0
        print *, "pr0 = ", pr0
        print *, "Tem0 = ", Tem0
        print *, "j0 = ", j0
        print *, "eta0 = ", eta0
        print *, "R0 = ", R0
        print *, "const1 = ", const1
        print *, "const2 = ", const2
        print *, "const3 = ", const3

        ! grid and index
        Nx=101
        Ny=101
        Nz=101
        xmin = -15.d0
        xmax = 15.d0
        ymin = -50.d0
        ymax = 50.d0
        zmin = -50.d0
        zmax = 50.d0
        
        call grid

        ! time and time step
        time = 0.0
        nstep = 0
        nmax = 10
        nout = 1    ! output frequency
        allocate(time_series(nmax))

        ! physical parameter
        allocate(ne(Nx, Ny, Nz), ni(Nx, Ny, Nz))
        allocate(vex(Nx, Ny, Nz), vey(Nx, Ny, Nz), vez(Nx, Ny, Nz), vix(Nx, Ny, Nz), viy(Nx, Ny, Nz), viz(Nx, Ny, Nz))
        allocate(Bsx(Nx, Ny, Nz), Bsy(Nx, Ny, Nz), Bsz(Nx, Ny, Nz), Bex(Nx, Ny, Nz), Bey(Nx, Ny, Nz), Bez(Nx, Ny, Nz))
        allocate(Esx(Nx, Ny, Nz), Esy(Nx, Ny, Nz), Esz(Nx, Ny, Nz), Eex(Nx, Ny, Nz), Eey(Nx, Ny, Nz), Eez(Nx, Ny, Nz))
        allocate(Bsx_pre(Nx,Ny,Nz), Bsy_pre(Nx,Ny,Nz), Bsz_pre(Nx,Ny,Nz))
        allocate(Esx_pre(Nx,Ny,Nz), Esy_pre(Nx,Ny,Nz), Esz_pre(Nx,Ny,Nz))
        allocate(Bx(Nx, Ny, Nz), By(Nx, Ny, Nz), Bz(Nx, Ny, Nz), Ex(Nx, Ny, Nz), Ey(Nx, Ny, Nz), Ez(Nx, Ny, Nz))
        allocate(jx(Nx, Ny, Nz), jy(Nx, Ny, Nz), jz(Nx, Ny, Nz))
        allocate(Te(Nx, Ny, Nz), Ti(Nx, Ny, Nz))
        allocate(pre(Nx, Ny, Nz), pri(Nx, Ny, Nz), eta(Nx, Ny, Nz))
        allocate(divB(Nx,Ny,Nz),divB_time(nmax),divE(Nx,Ny,Nz),divE_time(nmax))

        ne = 1.0d0
        ni = 1.0d0

        Bsx = 0.d0
        Bsy = 0.d0
        do i = 1, Nx
            Bsz(i,:,:) = Bh*tanh(x(i))
        end do
        Bex = 0.d0
        Bey = 0.d0
        Bez = 0.d0

        Esx = 0.d0
        Esy = 0.d0
        Esz = 0.d0
        Eex = 0.d0
        Eey = 0.d0
        Eez = 0.d0

        Bsx_pre = 0.d0
        Bsy_pre = 0.d0
        Bsz_pre = 0.d0
        Esx_pre = 0.d0
        Esy_pre = 0.d0
        Esz_pre = 0.d0

        call total_EMfield

        pre = B0**2/(2.0d0*mu0)*(Bh**2-Bsz**2)/(2*pr0)
        pri = B0**2/(2.0d0*mu0)*(Bh**2-Bsz**2)/(2*pr0)

        Te = pre/ne
        Ti = pri/ni

        vex = 0.d0
        vey = 0.0d0
        vez = 0.d0
        vix = 0.d0
        viy = 0.0d0
        viz = 0.d0
        
        eta = 0.0d0

        call current

        divB = 0.d0
        gamma = 5.0d0/3.0d0

        ! controling parameter
        is_abnormal_resistance = .false.
        
    end subroutine

    subroutine initialize_Harris_sheet
        implicit none

        ! parameters of Harris sheet -------------------------------
        real(kind=8) :: Lh = 1.0d0 ! Harris sheet half thickness
        real(kind=8) :: Bh = -1.0d0 ! Harris sheet magnetic field
        integer :: i
        ! --------------------------------------------------------

        ! constants and normalizing parameters
        x0 = 3.4d5
        t0 = 1.42344d-4
        n0 = 5.0d6
        v0 = x0 / t0
        E0 = m0 * v0 / (q0 * t0)
        B0 = m0 / (q0 * t0)
        pr0 = m0 * n0 * v0**2
        Tem0 = m0 * v0**2
        j0 = B0 / (mu0 * x0)
        eta0 = m0 / (q0**2 * n0 * t0)
        R0 = eta0 * q0 * n0 * j0
        const1 = mu0 * x0**2 * n0 * q0**2 / m0   ! used in calculating current density j
        const2 = t0**2 / (mu0 * eps0 * x0**2)    ! used in calculating electric field E
        const3 = R0 / (m0 * n0 * v0 / t0)   ! used in fraction term Re and Ri
        qe = -1.0d0
        qi = 1.0d0
        me = 1.0d0
        mi = 1864.0d0

        print *, "x0 = ", x0
        print *, "t0 = ", t0
        print *, "n0 = ", n0
        print *, "v0 = ", v0
        print *, "E0 = ", E0
        print *, "B0 = ", B0
        print *, "pr0 = ", pr0
        print *, "Tem0 = ", Tem0
        print *, "j0 = ", j0
        print *, "eta0 = ", eta0
        print *, "R0 = ", R0
        print *, "const1 = ", const1
        print *, "const2 = ", const2
        print *, "const3 = ", const3

        ! grid and index
        Nx=101
        Ny=101
        Nz=101
        xmin = -15.d0
        xmax = 15.d0
        ymin = -50.d0
        ymax = 50.d0
        zmin = -50.d0
        zmax = 50.d0

        call grid

        ! time and time step
        time = 0.0
        nstep = 0
        nmax = 10
        nout = 1
        allocate(time_series(nmax))

        ! physical parameter
        allocate(ne(Nx, Ny, Nz), ni(Nx, Ny, Nz))
        allocate(vex(Nx, Ny, Nz), vey(Nx, Ny, Nz), vez(Nx, Ny, Nz), vix(Nx, Ny, Nz), viy(Nx, Ny, Nz), viz(Nx, Ny, Nz))
        allocate(Bsx(Nx, Ny, Nz), Bsy(Nx, Ny, Nz), Bsz(Nx, Ny, Nz), Bex(Nx, Ny, Nz), Bey(Nx, Ny, Nz), Bez(Nx, Ny, Nz))
        allocate(Esx(Nx, Ny, Nz), Esy(Nx, Ny, Nz), Esz(Nx, Ny, Nz), Eex(Nx, Ny, Nz), Eey(Nx, Ny, Nz), Eez(Nx, Ny, Nz))
        allocate(Bsx_pre(Nx,Ny,Nz), Bsy_pre(Nx,Ny,Nz), Bsz_pre(Nx,Ny,Nz))
        allocate(Esx_pre(Nx,Ny,Nz), Esy_pre(Nx,Ny,Nz), Esz_pre(Nx,Ny,Nz))
        allocate(Bx(Nx, Ny, Nz), By(Nx, Ny, Nz), Bz(Nx, Ny, Nz), Ex(Nx, Ny, Nz), Ey(Nx, Ny, Nz), Ez(Nx, Ny, Nz))
        allocate(jx(Nx, Ny, Nz), jy(Nx, Ny, Nz), jz(Nx, Ny, Nz))
        allocate(Te(Nx, Ny, Nz), Ti(Nx, Ny, Nz))
        allocate(pre(Nx, Ny, Nz), pri(Nx, Ny, Nz), eta(Nx, Ny, Nz))
        allocate(divB(Nx,Ny,Nz),divB_time(nmax),divE(Nx,Ny,Nz),divE_time(nmax))

        do i = 1, Nx
            ne(i,:,:) = 1.0d0/cosh(x(i)/Lh)**2
            ni(i,:,:) = 1.0d0/cosh(x(i)/Lh)**2
        end do

        

        Bsx = 0.d0
        Bsy = 0.d0
        do i = 1, Nx
            Bsz(i,:,:) = Bh*tanh(x(i)/Lh)
        end do
        Bex = 0.d0
        Bey = 0.d0
        Bez = 0.d0
        Esx = 0.d0
        Esy = 0.d0
        Esz = 0.d0
        Eex = 0.d0
        Eey = 0.d0
        Eez = 0.d0

        Bsx_pre = 0.d0  ! used in leapfrog method
        Bsy_pre = 0.d0
        Bsz_pre = 0.d0
        Esx_pre = 0.d0
        Esy_pre = 0.d0
        Esz_pre = 0.d0

        call total_EMfield

        pre = B0**2/(2.0d0*mu0)*(Bh**2-Bsz**2)/(2*pr0)
        pri = B0**2/(2.0d0*mu0)*(Bh**2-Bsz**2)/(2*pr0)

        Te = pre/ne
        Ti = pri/ni

        vex = 0.d0
        vey = -2.0d0*Te/(qe*Bh*Lh)
        vez = 0.d0
        vix = 0.d0
        viy = -2.0d0*Ti/(qi*Bh*Lh)
        viz = 0.d0
        
        eta = 0.0d0

        call current

        divB = 0.d0
        gamma = 5.0d0/3.0d0

        ! controling parameter
        is_abnormal_resistance = .false.
        
    end subroutine

    subroutine stepon
        implicit none

        !call shear_flow

        call continuity_equation
        call momentum_equation
        call energy_equation
        call EMField_eqution
        
        call total_EMfield

        call abnormal_resistance
        call current
        call pressure

        call check
    end subroutine

    subroutine energy_equation
        implicit none
        real(kind=8), allocatable :: K1(:,:,:), K2(:,:,:), K3(:,:,:), K4(:,:,:), temp(:,:,:)
        real(kind=8), allocatable :: cdiffx(:,:,:), cdiffy(:,:,:), cdiffz(:,:,:)

        allocate(K1(Nx, Ny, Nz), K2(Nx, Ny, Nz), K3(Nx, Ny, Nz), K4(Nx, Ny, Nz), temp(Nx, Ny, Nz))
        allocate(cdiffx(Nx,Ny,Nz),cdiffy(Nx,Ny,Nz),cdiffz(Nx,Ny,Nz))

        ! Te first
        ! inner points use RK4
        K1 = -tau * (vex * central_difference_x(Te,hx) + vey * central_difference_y(Te,hy) + vez * central_difference_z(Te,hz))
        K2 = -tau * (vex * central_difference_x(Te+K1/2.0d0,hx) + vey * central_difference_y(Te+K1/2.0d0,hy) &
                     + vez * central_difference_z(Te+K1/2.0d0,hz))
        K3 = -tau * (vex * central_difference_x(Te+K2/2.0d0,hx) + vey * central_difference_y(Te+K2/2.0d0,hy) + & 
                        vez * central_difference_z(Te+K2/2.0d0,hz))
        K4 = -tau * (vex * central_difference_x(Te+K3,hx) + vey * central_difference_y(Te+K3,hy) & 
                        + vez * central_difference_z(Te+K3,hz))
        temp = Te + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(Te,hx)
        cdiffy = central_difference_y(Te,hy)
        cdiffz = central_difference_z(Te,hz)

        temp(2:4,:,:) = Te(2:4,:,:) - tau * (vex(2:4,:,:) * cdiffx(2:4,:,:) + vey(2:4,:,:) * cdiffy(2:4,:,:) & 
                                            + vez(2:4,:,:) * cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = Te(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:) * cdiffx(Nx-3:Nx-1,:,:) & 
                                        + vey(Nx-3:Nx-1,:,:) * cdiffy(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:) * cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = Te(:,2:4,:) - tau * (vex(:,2:4,:) * cdiffx(:,2:4,:) + vey(:,2:4,:) * cdiffy(:,2:4,:) & 
                                            + vez(:,2:4,:) * cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = Te(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:) * cdiffx(:,Ny-3:Ny-1,:) & 
                                        + vey(:,Ny-3:Ny-1,:) * cdiffy(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:) * cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = Te(:,:,2:4) - tau * (vex(:,:,2:4) * cdiffx(:,:,2:4) + vey(:,:,2:4) * cdiffy(:,:,2:4) & 
                                            + vez(:,:,2:4) * cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = Te(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1) * cdiffx(:,:,Nz-3:Nz-1) & 
                                        + vey(:,:,Nz-3:Nz-1) * cdiffy(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1) * cdiffz(:,:,Nz-3:Nz-1))

        ! add other terms
        temp = temp + tau * 2.0d0 / (3.0d0*ne) * (-pre * (central_difference_x(vex,hx) + central_difference_y(vey,hy) &
                                                            + central_difference_z(vez,hz)))

        ! boundary points
        call boundary(temp)

        ! update Te
        Te = temp

        ! Ti second
        ! inner points use RK4
        K1 = -tau * (vix * central_difference_x(Ti,hx) + viy * central_difference_y(Ti,hy) + viz * central_difference_z(Ti,hz))
        K2 = -tau * (vix * central_difference_x(Ti+K1/2.0d0,hx) + viy * central_difference_y(Ti+K1/2.0d0,hy) & 
                        + viz * central_difference_z(Ti+K1/2.0d0,hz))
        K3 = -tau * (vix * central_difference_x(Ti+K2/2.0d0,hx) + viy * central_difference_y(Ti+K2/2.0d0,hy) & 
                        + viz * central_difference_z(Ti+K2/2.0d0,hz))
        K4 = -tau * (vix * central_difference_x(Ti+K3,hx) + viy * central_difference_y(Ti+K3,hy) & 
                        + viz * central_difference_z(Ti+K3,hz))
        temp = Ti + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(Ti,hx)
        cdiffy = central_difference_y(Ti,hy)
        cdiffz = central_difference_z(Ti,hz)

        temp(2:4,:,:) = Ti(2:4,:,:) - tau * (vix(2:4,:,:) * cdiffx(2:4,:,:) + viy(2:4,:,:) * cdiffy(2:4,:,:) & 
                                                + viz(2:4,:,:) * cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = Ti(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:) * cdiffx(Nx-3:Nx-1,:,:) &
                                        + viy(Nx-3:Nx-1,:,:) * cdiffy(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:) * cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = Ti(:,2:4,:) - tau * (vix(:,2:4,:) * cdiffx(:,2:4,:) + viy(:,2:4,:) * cdiffy(:,2:4,:) & 
                                                + viz(:,2:4,:) * cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = Ti(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:) * cdiffx(:,Ny-3:Ny-1,:) & 
                                        + viy(:,Ny-3:Ny-1,:) * cdiffy(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:) * cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = Ti(:,:,2:4) - tau * (vix(:,:,2:4) * cdiffx(:,:,2:4) + viy(:,:,2:4) * cdiffy(:,:,2:4) &
                                                + viz(:,:,2:4) * cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = Ti(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1) * cdiffx(:,:,Nz-3:Nz-1) & 
                                        + viy(:,:,Nz-3:Nz-1) * cdiffy(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1) * cdiffz(:,:,Nz-3:Nz-1))

        ! add other terms
        temp = temp + tau * 2.0d0 / (3.0d0*ni) * (-pri * (central_difference_x(vix,hx) + central_difference_y(viy,hy) &
                                                            + central_difference_z(viz,hz)))

        ! boundary points
        call boundary(temp)

        ! update Ti
        Ti = temp

        deallocate(K1, K2, K3, K4, temp)
        deallocate(cdiffx, cdiffy, cdiffz)

    end subroutine

    subroutine continuity_equation
        implicit none
        real(kind=8), allocatable :: K1(:,:,:), K2(:,:,:), K3(:,:,:), K4(:,:,:), temp(:,:,:)
        real(kind=8), allocatable :: cdiffx(:,:,:), cdiffy(:,:,:), cdiffz(:,:,:)

        allocate(K1(Nx, Ny, Nz), K2(Nx, Ny, Nz), K3(Nx, Ny, Nz), K4(Nx, Ny, Nz), temp(Nx, Ny, Nz))
        allocate(cdiffx(Nx,Ny,Nz),cdiffy(Nx,Ny,Nz),cdiffz(Nx,Ny,Nz))

        ! ne first
        ! inner points use RK4
        K1 = -tau * (central_difference_x(ne*vex, hx)+central_difference_y(ne*vey, hy)+central_difference_z(ne*vez, hz))
        K2 = -tau * (central_difference_x((ne+K1/2.0d0)*vex, hx)+central_difference_y((ne+K1/2.0d0)*vey, hy) & 
                        +central_difference_z((ne+K1/2.0d0)*vez, hz))
        K3 = -tau * (central_difference_x((ne+K2/2.0d0)*vex, hx)+central_difference_y((ne+K2/2.0d0)*vey, hy) & 
                        +central_difference_z((ne+K2/2.0d0)*vez, hz))
        K4 = -tau * (central_difference_x((ne+K3)*vex, hx)+central_difference_y((ne+K3)*vey, hy) & 
                        +central_difference_z((ne+K3)*vez, hz))
        temp = ne + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(ne*vex,hx)
        cdiffy = central_difference_y(ne*vey,hy)
        cdiffz = central_difference_z(ne*vez,hz)

        temp(2:4,:,:) = ne(2:4,:,:) - tau * (cdiffx(2:4,:,:) + cdiffy(2:4,:,:) + cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = ne(Nx-3:Nx-1,:,:) - tau * (cdiffx(Nx-3:Nx-1,:,:) + cdiffy(Nx-3:Nx-1,:,:) + cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = ne(:,2:4,:) - tau * (cdiffx(:,2:4,:) + cdiffy(:,2:4,:) + cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = ne(:,Ny-3:Ny-1,:) - tau * (cdiffx(:,Ny-3:Ny-1,:) + cdiffy(:,Ny-3:Ny-1,:) + cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = ne(:,:,2:4) - tau * (cdiffx(:,:,2:4) + cdiffy(:,:,2:4) + cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = ne(:,:,Nz-3:Nz-1) - tau * (cdiffx(:,:,Nz-3:Nz-1) + cdiffy(:,:,Nz-3:Nz-1) + cdiffz(:,:,Nz-3:Nz-1))

        ! boundary points
        call boundary(temp)

        ! update ne
        ne = temp

        ! ni second
        ! inner points use RK4
        K1 = -tau * (central_difference_x(ni*vix, hx)+central_difference_y(ni*viy, hy)+central_difference_z(ni*viz, hz))
        K2 = -tau * (central_difference_x((ni+K1/2.0d0)*vix, hx)+central_difference_y((ni+K1/2.0d0)*viy, hy) & 
                        +central_difference_z((ni+K1/2.0d0)*viz, hz))
        K3 = -tau * (central_difference_x((ni+K2/2.0d0)*vix, hx)+central_difference_y((ni+K2/2.0d0)*viy, hy) & 
                        +central_difference_z((ni+K2/2.0d0)*viz, hz))
        K4 = -tau * (central_difference_x((ni+K3)*vix, hx)+central_difference_y((ni+K3)*viy, hy) & 
                        +central_difference_z((ni+K3)*viz, hz))
        temp = ni + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(ni*vix,hx)
        cdiffy = central_difference_y(ni*viy,hy)
        cdiffz = central_difference_z(ni*viz,hz)

        temp(2:4,:,:) = ni(2:4,:,:) - tau * (cdiffx(2:4,:,:) + cdiffy(2:4,:,:) + cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = ni(Nx-3:Nx-1,:,:) - tau * (cdiffx(Nx-3:Nx-1,:,:) + cdiffy(Nx-3:Nx-1,:,:) + cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = ni(:,2:4,:) - tau * (cdiffx(:,2:4,:) + cdiffy(:,2:4,:) + cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = ni(:,Ny-3:Ny-1,:) - tau * (cdiffx(:,Ny-3:Ny-1,:) + cdiffy(:,Ny-3:Ny-1,:) + cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = ni(:,:,2:4) - tau * (cdiffx(:,:,2:4) + cdiffy(:,:,2:4) + cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = ni(:,:,Nz-3:Nz-1) - tau * (cdiffx(:,:,Nz-3:Nz-1) + cdiffy(:,:,Nz-3:Nz-1) + cdiffz(:,:,Nz-3:Nz-1))
        ! boundary points
        call boundary(temp)

        ! update ni
        ni = temp

        deallocate(K1, K2, K3, K4, temp)
        deallocate(cdiffx, cdiffy, cdiffz)
    end subroutine

    subroutine momentum_equation
        implicit none
        real(kind=8), allocatable :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
        real(kind=8), allocatable :: K1(:,:,:), K2(:,:,:), K3(:,:,:), K4(:,:,:)
        ! integer :: i, j, k
        ! real(kind=8) :: cdiffz(Nx,Ny,Nz)

        allocate(temp1(Nx,Ny,Nz), temp2(Nx,Ny,Nz), temp3(Nx,Ny,Nz))
        allocate(K1(Nx,Ny,Nz), K2(Nx,Ny,Nz), K3(Nx,Ny,Nz), K4(Nx,Ny,Nz))

        call fraction_term

        ! electron
        ! inner points use RK4
        ! vex
        K1 = -tau * (vex * central_difference_x(vex, hx) + vey * central_difference_y(vex, hy) & 
                        + vez * central_difference_z(vex, hz))
        K2 = -tau * ((vex+K1/2.0d0) * central_difference_x(vex+K1/2.0d0, hx) + vey * central_difference_y(vex+K1/2.0d0, hy) & 
                        + vez * central_difference_z(vex+K1/2.0d0, hz))
        K3 = -tau * ((vex+K2/2.0d0) * central_difference_x(vex+K2/2.0d0, hx) + vey * central_difference_y(vex+K2/2.0d0, hy) & 
                        + vez * central_difference_z(vex+K2/2.0d0, hz))
        K4 = -tau * ((vex+K3) * central_difference_x(vex+K3, hx) + vey * central_difference_y(vex+K3, hy) & 
                        + vez * central_difference_z(vex+K3, hz))
        temp1 = vex + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! vey
        K1 = -tau * (vex * central_difference_x(vey, hx) + vey * central_difference_y(vey, hy) & 
                        + vez * central_difference_z(vey, hz))
        K2 = -tau * (vex * central_difference_x(vey+K1/2.0d0, hx) + (vey+K1/2.0d0) * central_difference_y(vey+K1/2.0d0, hy) & 
                        + vez * central_difference_z(vey+K1/2.0d0, hz))
        K3 = -tau * (vex * central_difference_x(vey+K2/2.0d0, hx) + (vey+K2/2.0d0) * central_difference_y(vey+K2/2.0d0, hy) & 
                        + vez * central_difference_z(vey+K2/2.0d0, hz))
        K4 = -tau * (vex * central_difference_x(vey+K3, hx) + (vey+K3) * central_difference_y(vey+K3, hy) & 
                        + vez * central_difference_z(vey+K3, hz))
        temp2 = vey + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! vez
        K1 = -tau * (vex * central_difference_x(vez, hx) + vey * central_difference_y(vez, hy) & 
                        + vez * central_difference_z(vez, hz))
        K2 = -tau * (vex * central_difference_x(vez+K1/2.0d0, hx) + vey * central_difference_y(vez+K1/2.0d0, hy) & 
                        + (vez+K1/2.0d0) * central_difference_z(vez+K1/2.0d0, hz))
        K3 = -tau * (vex * central_difference_x(vez+K2/2.0d0, hx) + vey * central_difference_y(vez+K2/2.0d0, hy) & 
                        + (vez+K2/2.0d0) * central_difference_z(vez+K2/2.0d0, hz))
        K4 = -tau * (vex * central_difference_x(vez+K3, hx) + vey * central_difference_y(vez+K3, hy) & 
                        + (vez+K3) * central_difference_z(vez+K3, hz))
        temp3 = vez + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        

        ! subouter points use Euler
        call momentum_subouter(temp1, temp2, temp3, 1)

        ! add other terms
        temp1 = temp1 + tau / (me * ne) * (ne * qe * (Ex + (vey * Bz - vez * By)) - central_difference_x(pre, hx) + const3 * Rex)
        temp2 = temp2 + tau / (me * ne) * (ne * qe * (Ey + (vez * Bx - vex * Bz)) - central_difference_y(pre, hy) + const3 * Rey)
        temp3 = temp3 + tau / (me * ne) * (ne * qe * (Ez + (vex * By - vey * Bx)) - central_difference_z(pre, hz) + const3 * Rez)

        ! boundary points
        call boundary(temp1)
        call boundary(temp2)
        call boundary(temp3)

        ! update
        vex = temp1
        vey = temp2
        vez = temp3

        ! ion
        ! inner points use RK4
        ! vix
        K1 = -tau * (vix * central_difference_x(vix, hx) + viy * central_difference_y(vix, hy) & 
                        + viz * central_difference_z(vix, hz))
        K2 = -tau * ((vix+K1/2.0d0) * central_difference_x(vix+K1/2.0d0, hx) + viy * central_difference_y(vix+K1/2.0d0, hy) & 
                        + viz * central_difference_z(vix+K1/2.0d0, hz))
        K3 = -tau * ((vix+K2/2.0d0) * central_difference_x(vix+K2/2.0d0, hx) + viy * central_difference_y(vix+K2/2.0d0, hy) & 
                        + viz * central_difference_z(vix+K2/2.0d0, hz))
        K4 = -tau * ((vix+K3) * central_difference_x(vix+K3, hx) + viy * central_difference_y(vix+K3, hy) & 
                        + viz * central_difference_z(vix+K3, hz))
        temp1 = vix + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! viy
        K1 = -tau * (vix * central_difference_x(viy, hx) + viy * central_difference_y(viy, hy) & 
                        + viz * central_difference_z(viy, hz))
        K2 = -tau * (vix * central_difference_x(viy+K1/2.0d0, hx) + (viy+K1/2.0d0) * central_difference_y(viy+K1/2.0d0, hy) & 
                        + viz * central_difference_z(viy+K1/2.0d0, hz))
        K3 = -tau * (vix * central_difference_x(viy+K2/2.0d0, hx) + (viy+K2/2.0d0) * central_difference_y(viy+K2/2.0d0, hy) & 
                        + viz * central_difference_z(viy+K2/2.0d0, hz))
        K4 = -tau * (vix * central_difference_x(viy+K3, hx) + (viy+K3) * central_difference_y(viy+K3, hy) & 
                        + viz * central_difference_z(viy+K3, hz))
        temp2 = viy + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! viz
        K1 = -tau * (vix * central_difference_x(viz, hx) + viy * central_difference_y(viz, hy) & 
                        + viz * central_difference_z(viz, hz))
        K2 = -tau * (vix * central_difference_x(viz+K1/2.0d0, hx) + viy * central_difference_y(viz+K1/2.0d0, hy) & 
                        + (viz+K1/2.0d0) * central_difference_z(viz+K1/2.0d0, hz))
        K3 = -tau * (vix * central_difference_x(viz+K2/2.0d0, hx) + viy * central_difference_y(viz+K2/2.0d0, hy) & 
                        + (viz+K2/2.0d0) * central_difference_z(viz+K2/2.0d0, hz))
        K4 = -tau * (vix * central_difference_x(viz+K3, hx) + viy * central_difference_y(viz+K3, hy) & 
                        + (viz+K3) * central_difference_z(viz+K3, hz))
        temp3 = viz + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! open(unit=104,file="temp31.dat",status='unknown',form='formatted')
        ! write(104,*)'TITLE=Debug'
        ! write(104,*) 'VARIABLES="x" "y" "z" "data"'
        ! write(104,"('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')") Nx,Ny,Nz
        ! write(104,"(4(1x,e11.4))") (((x(i), y(j), z(k), temp3(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        ! close(104)
        

        ! subouter points use Euler
        call momentum_subouter(temp1, temp2, temp3, 2)
        
        ! open(unit=104,file="temp32.dat",status='unknown',form='formatted')
        ! write(104,*)'TITLE=Debug'
        ! write(104,*) 'VARIABLES="x" "y" "z" "data"'
        ! write(104,"('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')") Nx,Ny,Nz
        ! write(104,"(4(1x,e11.4))") (((x(i), y(j), z(k), temp3(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        ! close(104)

        ! add other terms
        temp1 = temp1 + tau / (mi * ni) * (ni * qi * (Ex + (viy * Bz - viz * By)) - central_difference_x(pri, hx) + const3 * Rix)
        temp2 = temp2 + tau / (mi * ni) * (ni * qi * (Ey + (viz * Bx - vix * Bz)) - central_difference_y(pri, hy) + const3 * Riy)
        temp3 = temp3 + tau / (mi * ni) * (ni * qi * (Ez + (vix * By - viy * Bx)) - central_difference_z(pri, hz) + const3 * Riz)

        ! open(unit=104,file="temp33.dat",status='unknown',form='formatted')
        ! write(104,*)'TITLE=Debug'
        ! write(104,*) 'VARIABLES="x" "y" "z" "data"'
        ! write(104,"('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')") Nx,Ny,Nz
        ! write(104,"(4(1x,e11.4))") (((x(i), y(j), z(k), temp3(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        ! close(104)

        ! cdiffz = central_difference_z(pri, hz)

        ! open(unit=104,file="detail.dat",status='unknown',form='formatted')
        ! write(104,*)'TITLE=Debug'
        ! write(104,*) 'VARIABLES="x" "y" "z" "vixBy" "viyBx" "cdiffz" "pri"'
        ! write(104,"('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')") Nx,Ny,Nz
        ! write(104,"(7(1x,e11.4))") & 
        !     (((x(i), y(j), z(k), vix(i,j,k)*By(i,j,k), viy(i,j,k)*Bx(i,j,k), cdiffz(i,j,k), pri(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        ! close(104)

        ! boundary points
        call boundary(temp1)
        call boundary(temp2)
        call boundary(temp3)

        ! open(unit=104,file="temp34.dat",status='unknown',form='formatted')
        ! write(104,*)'TITLE=Debug'
        ! write(104,*) 'VARIABLES="x" "y" "z" "data"'
        ! write(104,"('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')") Nx,Ny,Nz
        ! write(104,"(4(1x,e11.4))") (((x(i), y(j), z(k), temp3(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        ! close(104)

        ! update
        vix = temp1
        viy = temp2
        viz = temp3

        deallocate(temp1, temp2, temp3)
        deallocate(K1, K2, K3, K4)
    end subroutine

    subroutine fraction_term
        implicit none

        Rex = eta * (-qe) * ne * jx
        Rey = eta * (-qe) * ne * jy
        Rez = eta * (-qe) * ne * jz

        Rix = -Rex
        Riy = -Rey
        Riz = -Rez
    end subroutine

    subroutine momentum_subouter(vx, vy, vz, species)
        implicit none
        real(kind=8) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        integer :: species  ! 1 for electron, 2 for ion
        real(kind=8), allocatable :: cdiffxx(:,:,:), cdiffxy(:,:,:), cdiffxz(:,:,:)
        real(kind=8), allocatable :: cdiffyx(:,:,:), cdiffyy(:,:,:), cdiffyz(:,:,:)
        real(kind=8), allocatable :: cdiffzx(:,:,:), cdiffzy(:,:,:), cdiffzz(:,:,:)

        allocate(cdiffxx(Nx,Ny,Nz), cdiffxy(Nx,Ny,Nz), cdiffxz(Nx,Ny,Nz))
        allocate(cdiffyx(Nx,Ny,Nz), cdiffyy(Nx,Ny,Nz), cdiffyz(Nx,Ny,Nz))
        allocate(cdiffzx(Nx,Ny,Nz), cdiffzy(Nx,Ny,Nz), cdiffzz(Nx,Ny,Nz))

        if (species==1) then    ! electron
            cdiffxx = central_difference_x(vex, hx)
            cdiffxy = central_difference_x(vey, hx)
            cdiffxz = central_difference_x(vez, hx)
            cdiffyx = central_difference_y(vex, hy)
            cdiffyy = central_difference_y(vey, hy)
            cdiffyz = central_difference_y(vez, hy)
            cdiffzx = central_difference_z(vex, hz)
            cdiffzy = central_difference_z(vey, hz)
            cdiffzz = central_difference_z(vez, hz)

            vx(2:4,:,:) = vex(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxx(2:4,:,:)+ vey(2:4,:,:)*cdiffyx(2:4,:,:) & 
                                                    + vez(2:4,:,:)*cdiffzx(2:4,:,:))
            vx(Nx-3:Nx-1,:,:) = vex(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxx(Nx-3:Nx-1,:,:) & 
                                            + vey(Nx-3:Nx-1,:,:)*cdiffyx(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzx(Nx-3:Nx-1,:,:))
            vx(:,2:4,:) = vex(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxx(:,2:4,:)+ vey(:,2:4,:)*cdiffyx(:,2:4,:) & 
                                                    + vez(:,2:4,:)*cdiffzx(:,2:4,:))
            vx(:,Ny-3:Ny-1,:) = vex(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxx(:,Ny-3:Ny-1,:) & 
                                            + vey(:,Ny-3:Ny-1,:)*cdiffyx(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzx(:,Ny-3:Ny-1,:))
            vx(:,:,2:4) = vex(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxx(:,:,2:4)+ vey(:,:,2:4)*cdiffyx(:,:,2:4) & 
                                                    + vez(:,:,2:4)*cdiffzx(:,:,2:4))
            vx(:,:,Nz-3:Nz-1) = vex(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxx(:,:,Nz-3:Nz-1) & 
                                            + vey(:,:,Nz-3:Nz-1)*cdiffyx(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzx(:,:,Nz-3:Nz-1))
            
            vy(2:4,:,:) = vey(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxy(2:4,:,:)+ vey(2:4,:,:)*cdiffyy(2:4,:,:) & 
                                                    + vez(2:4,:,:)*cdiffzy(2:4,:,:))
            vy(Nx-3:Nx-1,:,:) = vey(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxy(Nx-3:Nx-1,:,:) & 
                                            + vey(Nx-3:Nx-1,:,:)*cdiffyy(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzy(Nx-3:Nx-1,:,:))
            vy(:,2:4,:) = vey(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxy(:,2:4,:)+ vey(:,2:4,:)*cdiffyy(:,2:4,:) & 
                                                    + vez(:,2:4,:)*cdiffzy(:,2:4,:))
            vy(:,Ny-3:Ny-1,:) = vey(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxy(:,Ny-3:Ny-1,:) & 
                                            + vey(:,Ny-3:Ny-1,:)*cdiffyy(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzy(:,Ny-3:Ny-1,:))
            vy(:,:,2:4) = vey(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxy(:,:,2:4)+ vey(:,:,2:4)*cdiffyy(:,:,2:4) & 
                                                    + vez(:,:,2:4)*cdiffzy(:,:,2:4))
            vy(:,:,Nz-3:Nz-1) = vey(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxy(:,:,Nz-3:Nz-1) & 
                                            + vey(:,:,Nz-3:Nz-1)*cdiffyy(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzy(:,:,Nz-3:Nz-1))
            
            vz(2:4,:,:) = vez(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxz(2:4,:,:)+ vey(2:4,:,:)*cdiffyz(2:4,:,:) & 
                                                    + vez(2:4,:,:)*cdiffzz(2:4,:,:))
            vz(Nx-3:Nx-1,:,:) = vez(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxz(Nx-3:Nx-1,:,:) & 
                                            + vey(Nx-3:Nx-1,:,:)*cdiffyz(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzz(Nx-3:Nx-1,:,:))
            vz(:,2:4,:) = vez(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxz(:,2:4,:)+ vey(:,2:4,:)*cdiffyz(:,2:4,:) & 
                                                    + vez(:,2:4,:)*cdiffzz(:,2:4,:))
            vz(:,Ny-3:Ny-1,:) = vez(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxz(:,Ny-3:Ny-1,:) & 
                                            + vey(:,Ny-3:Ny-1,:)*cdiffyz(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzz(:,Ny-3:Ny-1,:))
            vz(:,:,2:4) = vez(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxz(:,:,2:4)+ vey(:,:,2:4)*cdiffyz(:,:,2:4) &
                                                    + vez(:,:,2:4)*cdiffzz(:,:,2:4))
            vz(:,:,Nz-3:Nz-1) = vez(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxz(:,:,Nz-3:Nz-1) & 
                                            + vey(:,:,Nz-3:Nz-1)*cdiffyz(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzz(:,:,Nz-3:Nz-1))
   
        else if (species==2) then   ! ion
            cdiffxx = central_difference_x(vix, hx)
            cdiffxy = central_difference_x(viy, hx)
            cdiffxz = central_difference_x(viz, hx)
            cdiffyx = central_difference_y(vix, hy)
            cdiffyy = central_difference_y(viy, hy)
            cdiffyz = central_difference_y(viz, hy)
            cdiffzx = central_difference_z(vix, hz)
            cdiffzy = central_difference_z(viy, hz)
            cdiffzz = central_difference_z(viz, hz)

            vx(2:4,:,:) = vix(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxx(2:4,:,:)+ viy(2:4,:,:)*cdiffyx(2:4,:,:) & 
                                                    + viz(2:4,:,:)*cdiffzx(2:4,:,:))
            vx(Nx-3:Nx-1,:,:) = vix(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxx(Nx-3:Nx-1,:,:) & 
                                            + viy(Nx-3:Nx-1,:,:)*cdiffyx(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzx(Nx-3:Nx-1,:,:))
            vx(:,2:4,:) = vix(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxx(:,2:4,:)+ viy(:,2:4,:)*cdiffyx(:,2:4,:) & 
                                                    + viz(:,2:4,:)*cdiffzx(:,2:4,:))
            vx(:,Ny-3:Ny-1,:) = vix(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxx(:,Ny-3:Ny-1,:) & 
                                            + viy(:,Ny-3:Ny-1,:)*cdiffyx(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzx(:,Ny-3:Ny-1,:))
            vx(:,:,2:4) = vix(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxx(:,:,2:4)+ viy(:,:,2:4)*cdiffyx(:,:,2:4) & 
                                                    + viz(:,:,2:4)*cdiffzx(:,:,2:4))
            vx(:,:,Nz-3:Nz-1) = vix(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxx(:,:,Nz-3:Nz-1) & 
                                            + viy(:,:,Nz-3:Nz-1)*cdiffyx(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzx(:,:,Nz-3:Nz-1))
            
            vy(2:4,:,:) = viy(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxy(2:4,:,:)+ viy(2:4,:,:)*cdiffyy(2:4,:,:) & 
                                                    + viz(2:4,:,:)*cdiffzy(2:4,:,:))
            vy(Nx-3:Nx-1,:,:) = viy(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxy(Nx-3:Nx-1,:,:) & 
                                            + viy(Nx-3:Nx-1,:,:)*cdiffyy(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzy(Nx-3:Nx-1,:,:))
            vy(:,2:4,:) = viy(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxy(:,2:4,:)+ viy(:,2:4,:)*cdiffyy(:,2:4,:) & 
                                                    + viz(:,2:4,:)*cdiffzy(:,2:4,:))
            vy(:,Ny-3:Ny-1,:) = viy(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxy(:,Ny-3:Ny-1,:) & 
                                            + viy(:,Ny-3:Ny-1,:)*cdiffyy(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzy(:,Ny-3:Ny-1,:))
            vy(:,:,2:4) = viy(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxy(:,:,2:4)+ viy(:,:,2:4)*cdiffyy(:,:,2:4) & 
                                                    + viz(:,:,2:4)*cdiffzy(:,:,2:4))
            vy(:,:,Nz-3:Nz-1) = viy(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxy(:,:,Nz-3:Nz-1) & 
                                            + viy(:,:,Nz-3:Nz-1)*cdiffyy(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzy(:,:,Nz-3:Nz-1))
            
            vz(2:4,:,:) = viz(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxz(2:4,:,:)+ viy(2:4,:,:)*cdiffyz(2:4,:,:) & 
                                                    + viz(2:4,:,:)*cdiffzz(2:4,:,:))
            vz(Nx-3:Nx-1,:,:) = viz(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxz(Nx-3:Nx-1,:,:) & 
                                            + viy(Nx-3:Nx-1,:,:)*cdiffyz(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzz(Nx-3:Nx-1,:,:))
            vz(:,2:4,:) = viz(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxz(:,2:4,:)+ viy(:,2:4,:)*cdiffyz(:,2:4,:) & 
                                                    + viz(:,2:4,:)*cdiffzz(:,2:4,:))
            vz(:,Ny-3:Ny-1,:) = viz(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxz(:,Ny-3:Ny-1,:) & 
                                            + viy(:,Ny-3:Ny-1,:)*cdiffyz(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzz(:,Ny-3:Ny-1,:))
            vz(:,:,2:4) = viz(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxz(:,:,2:4)+ viy(:,:,2:4)*cdiffyz(:,:,2:4) &
                                            + viz(:,:,2:4)*cdiffzz(:,:,2:4))
            vz(:,:,Nz-3:Nz-1) = viz(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxz(:,:,Nz-3:Nz-1) & 
                                            + viy(:,:,Nz-3:Nz-1)*cdiffyz(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzz(:,:,Nz-3:Nz-1))

        end if

        deallocate(cdiffxx, cdiffxy, cdiffxz, cdiffyx, cdiffyy, cdiffyz, cdiffzx, cdiffzy, cdiffzz)
    end subroutine

    subroutine EMField_eqution
        implicit none

        real(kind=8), allocatable :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)   ! for Magnetic field
        real(kind=8), allocatable :: temp4(:,:,:), temp5(:,:,:), temp6(:,:,:)   ! for Electric field
        
        allocate(temp1(Nx,Ny,Nz), temp2(Nx,Ny,Nz), temp3(Nx,Ny,Nz))
        allocate(temp4(Nx,Ny,Nz), temp5(Nx,Ny,Nz), temp6(Nx,Ny,Nz))

        if (nstep==1) then
            ! Euler scheme
            temp1 = Bsx - tau * (central_difference_y(Esz,hy) - central_difference_z(Esy,hz))
            temp2 = Bsy - tau * (central_difference_z(Esx,hz) - central_difference_x(Esz,hx))
            temp3 = Bsz - tau * (central_difference_x(Esy,hx) - central_difference_y(Esx,hy))

            temp4 = Esx + tau * const2 * (central_difference_y(Bsz,hy) - central_difference_z(Bsy,hz) - Jx)
            temp5 = Esy + tau * const2 * (central_difference_z(Bsx,hz) - central_difference_x(Bsz,hx) - Jy)
            temp6 = Esz + tau * const2 * (central_difference_x(Bsy,hx) - central_difference_y(Bsx,hy) - Jz)
            
            !call record_debug(central_difference_x(Bsy,hx), "central_difference_x(Bsy,hx).dat")

            !call record_debug(central_difference_y(Bsx,hy), "central_difference_y(Bsx,hy).dat")

            !call record_debug(Jz, "Jz.dat")

            !call record_debug(temp6, "temp6.dat")

        else
            ! leapfrog scheme
            temp1 = Bsx_pre - 2.0d0 * tau * (central_difference_y(Esz,hy) - central_difference_z(Esy,hz))
            temp2 = Bsy_pre - 2.0d0 * tau * (central_difference_z(Esx,hz) - central_difference_x(Esz,hx))
            temp3 = Bsz_pre - 2.0d0 * tau * (central_difference_x(Esy,hx) - central_difference_y(Esx,hy))

            temp4 = Esx_pre + 2.0d0 * tau * const2 * (central_difference_y(Bsz,hy) - central_difference_z(Bsy,hz) - Jx)
            temp5 = Esy_pre + 2.0d0 * tau * const2 * (central_difference_z(Bsx,hz) - central_difference_x(Bsz,hx) - Jy)
            temp6 = Esz_pre + 2.0d0 * tau * const2 * (central_difference_x(Bsy,hx) - central_difference_y(Bsx,hy) - Jz)
        end if

        call boundary(temp1)
        call boundary(temp2)
        call boundary(temp3)
        call boundary(temp4)
        call boundary(temp5)
        call boundary(temp6)

        ! update
        Bsx_pre = Bsx
        Bsy_pre = Bsy
        Bsz_pre = Bsz

        Esx_pre = Esx
        Esy_pre = Esy
        Esz_pre = Esz

        Bsx = temp1
        Bsy = temp2
        Bsz = temp3

        Esx = temp4
        Esy = temp5
        Esz = temp6

        deallocate(temp1, temp2, temp3)
        deallocate(temp4, temp5, temp6)

    end subroutine

    subroutine record_debug(a, output)
        implicit none

        real(kind=8) :: a(:,:,:)

        integer :: i, j, k
        character(len=32) :: output
        character(len=50) :: format

        open(unit=103,file=output,status='unknown',form='formatted')
        ! write(103,*)'TITLE=Debug'
        ! write(103,*) 'VARIABLES="x" "y" "data"'
        ! format1 = "('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')"
        ! write(103,format1) Nx,Ny,Nz
        format = "(4(1x,e11.4))"
        write(103,format) (((x(i), y(j), z(k), a(i,j,k), i=1, Nx), j=1, Ny), k=1, Nz)
        close(103)

    end subroutine

    subroutine total_EMfield
        implicit none

        Bx = Bsx + Bex
        By = Bsy + Bey
        Bz = Bsz + Bez
        Ex = Esx + Eex
        Ey = Esy + Eey
        Ez = Esz + Eez

    end subroutine

    subroutine boundary(a)
        ! use linear extrapolation
        implicit none
        real(kind=8), intent(inout) :: a(:,:,:)
        integer :: size1, size2, size3

        size1 = size(a, 1)
        size2 = size(a, 2)
        size3 = size(a, 3)

        a(1,:,:) = 2.0d0 * a(2,:,:) - a(3,:,:)
        a(size1,:,:) = 2.0d0 * a(size1-1,:,:) - a(size1-2,:,:)
        a(:,1,:) = 2.0d0 * a(:,2,:) - a(:,3,:)
        a(:,size2,:) = 2.0d0 * a(:,size2-1,:) - a(:,size2-2,:)
        a(:,:,1) = 2.0d0 * a(:,:,2) - a(:,:,3)
        a(:,:,size3) = 2.0d0 * a(:,:,size3-1) - a(:,:,size3-2)
    end subroutine

    subroutine check
        ! check the divergence of B
        implicit none

        divB = central_difference_x(Bsx,hx) + central_difference_y(Bsy,hy) + central_difference_z(Bsz,hz)

        divB_max = maxval(maxval(maxval(divB,3),2),1)

        divE = central_difference_x(Esx,hx) + central_difference_y(Esy,hy) + central_difference_z(Esz,hz)

        divE_max = maxval(maxval(maxval(divE,3),2),1)

        divB_time(nstep) = divB_max
        divE_time(nstep) = divE_max

        print *, "divB= ", divB_max
        print *, "divE= ", divE_max
    end subroutine

    subroutine setdt

        implicit none

        real(kind=8), allocatable ::  temp1(:,:,:), temp2(:,:,:)
        real(kind=8) :: dxyz, dt1, dt2

        allocate(temp1(Nx,Ny,Nz), temp2(Nx,Ny,Nz))
        !call foreta(time, 1)
        call pressure

        dxyz = .5*hx*hy*hz/sqrt((hx**2*hy**2 + hy**2*hz**2 + hx**2*hz**2))

        temp1 = dxyz / (sqrt(vex**2 + vey**2 + vez**2) & 
                        + sqrt((B0**2/(2*mu0)*(Bx**2 + By**2 + Bz**2) + gamma*pr0*pre)/(m0*n0*me*ne))/v0)
        temp2 = dxyz / (sqrt(vix**2 + viy**2 + viz**2) & 
                        + sqrt((B0**2/(2*mu0)*(Bx**2 + By**2 + Bz**2) + gamma*pr0*pri)/(m0*n0*mi*ni))/v0)

        dt1 = minval(minval(minval(temp1,3),2),1)
        dt2 = minval(minval(minval(temp2,3),2),1)

        !tau = 0.5*min(dt1, dt2)
        tau = 0.01d0
        
        deallocate(temp1, temp2)
    end subroutine

    subroutine current
        ! calculate the current
        implicit none

        jx = const1 * (ne*qe*vex + ni*qi*vix)
        jy = const1 * (ne*qe*vey + ni*qi*viy)
        jz = const1 * (ne*qe*vez + ni*qi*viz)
        
    end subroutine



    subroutine grid

        implicit none

        integer :: i, j, k

        allocate(x(Nx), y(Ny), z(Nz))

        hx = (xmax - xmin)/(Nx-1)
        hy = (ymax - ymin)/(Ny-1)
        hz = (zmax - zmin)/(Nz-1)

        do i = 1, Nx
            x(i) = xmin + (i - 1)*hx
        end do

        do j = 1, Ny
            y(j) = ymin + (j - 1)*hy
        end do

        do k = 1, Nz
            z(k) = zmin + (k - 1)*hz
        end do

    end subroutine

    subroutine smooth(a, weight)

        implicit none
        real(kind=8) ::  a(:,:,:)
        real(kind=8) :: weight   ! weight of the central element

        real(kind=8), allocatable :: average(:,:,:)

        allocate(average(Nx,Ny,Nz))

        average(2:Nx-1,2:Ny-1,2:Nz-1) = (a(3:Nx,:,:) + a(1:Nx-2,:,:) + a(:,3:Ny,:) + a(:,1:Ny-2,:) &
                                         + a(:,:,3:Nz) + a(:,:,1:Nz-2)) / 6.0d0
        a = weight*a + (1 - weight)*average

    
        call boundary(a)

        deallocate(average)

    end subroutine

    subroutine pressure
        ! calculate the pressure of electron and ion
        implicit none

        pre = ne * Te
        pri = ni * Ti
    end subroutine

    

    subroutine abnormal_resistance
        implicit none

        if(.not. is_abnormal_resistance) return

    end subroutine

    
    subroutine shear_flow
        implicit none

        real(kind=8) :: vshear, tshear_start, tshear_end
        real(kind=8) :: xs, zs, lsx, lsz, tao_s
        integer :: i, k

        vshear = 3.2d-3  ! the velocity of the shear flow
        tshear_start = 0.0d0
        tshear_end = 10.0d0

        xs = 0.0d0  ! the central position of the shear flow
        zs = 0.0d0  ! the central position of the shear flow
        lsx = 1.0d0 ! the characteristic length of the shear flow
        lsz = 1.0d0 ! the characteristic length of the shear flow
        tao_s = 0.1d0  ! the characteristic time of the shear flow

        if (time>tshear_start .and. time<tshear_end) then
            do k = 1, Nz
                do i = 1, Nx
                    vey(i,:,k) = vey(i,:,k) - vshear * exp(-(((x(i) - xs)/lsx)**2  + ((z(k) - zs)/lsz)**2)) &
                                                 * (tanh((time+tau-tshear_start)/tao_s)-tanh((time-tshear_start)/tao_s))
                    viy(i,:,k) = viy(i,:,k) + vshear * exp(-(((x(i) - xs)/lsx)**2  + ((z(k) - zs)/lsz)**2)) &
                                                 * (tanh((time+tau-tshear_start)/tao_s)-tanh((time-tshear_start)/tao_s))
                end do
            end do
        end if

    end subroutine

    subroutine record
        ! record the data
        implicit none

        integer :: i, j, k
        character(len=50) :: output, format1, format2

        output = 'data_set_' // int2str(nstep) // '.dat'

        open(unit=101,file=output,status='unknown',form='formatted')
        write(101,*) 'TITLE=Two_Fluid_MHD'
        write(101,*) 'VARIABLES="x" "y" "z" "ne" "ni" "vex" "vey" "vez" "vix" "viy" "viz" &
                                &"Bx" "By" "Bz" "Ex" "Ey" "Ez" "Te" "Ti" "pe" "pi" "jx" "jy" "jz" '
        format1 = "('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')"
        write(101,format1) Nx,Ny,Nz

        format2 = "(24(1x,e11.4))"
        write(101,format2) (((x(i),y(j),z(k),ne(i,j,k),ni(i,j,k),vex(i,j,k),vey(i,j,k),vez(i,j,k), &
                        vix(i,j,k),viy(i,j,k),viz(i,j,k),Bx(i,j,k), By(i,j,k), Bz(i,j,k), &
                        Ex(i,j,k), Ey(i,j,k), Ez(i,j,k), Te(i,j,k), Ti(i,j,k), pre(i,j,k), pri(i,j,k), &
                        jx(i,j,k), jy(i,j,k), jz(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
    end subroutine

    subroutine record_divBE
        implicit none

        integer :: i
        character(len=50) :: output, format

        output = 'divBE' // '.dat'

        open(unit=102,file=output,status='unknown',form='formatted')
        write(102,*)'TITLE=Divergence_of_B_and_E'
        write(102,*) 'VARIABLES="divB" "divE" "time"'
        format = "(3(1x,e11.4))"
        write(102,format) (divB_time(i), divE_time(i), time_series(i), i=1, nstep)
    end subroutine

    subroutine destroy
        ! end of programm
        implicit none

        deallocate(x, y, z)
        deallocate(ne, ni)
        deallocate(vex, vey, vez, vix, viy, viz)
        deallocate(Bsx, Bsy, Bsz, Bex, Bey, Bez)
        deallocate(Esx, Esy, Esz, Eex, Eey, Eez)
        deallocate(Bsx_pre, Bsy_pre, Bsz_pre)
        deallocate(Esx_pre, Esy_pre, Esz_pre)
        deallocate(Bx, By, Bz, Ex, Ey, Ez)
        deallocate(jx, jy, jz)
        deallocate(Te, Ti)
        deallocate(pre, pri, eta)
        deallocate(divB)
    end subroutine

end module Custom_subroutines
