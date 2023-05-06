module custom_subroutines
    use define_variables
    use custom_functions
    implicit none

contains

    subroutine initialize
        implicit none
        ! constants and normalizing parameters
        v0 = x0 / t0
        E0 = m0 * v0 / (q0 * t0)
        B0 = m0 / (q0 * t0)
        pr0 = m0 * n0 * v0**2
        Tem0 = m0 * v0**2
        j0 = B0 / (mu0 * x0)
        const1 = mu0 * x0**2 * n0 * q0**2 / m0   ! used in calculating current density j
        const2 = t0**2 / (mu0 * eps0 * x0**2)    ! used in calculating electric field E

        
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
        allocate(x(Nx), y(Ny), z(Nz))
        call grid

        ! time and time step
        time = 0.0
        nstep = 0

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
        allocate(divB(Nx,Ny,Nz))
        ne = 1.d0
        ni = 1.d0
        vex = 0.d0
        vey = 0.d0
        vez = 0.d0
        vix = 0.d0
        viy = 0.d0
        viz = 0.d0
        Bsx = 0.d0
        Bsy = 0.d0
        Bsz = 0.d0
        Bex = 0.d0
        Bey = 0.d0
        Bez = 0.d0
        Esx = 0.d0
        Esy = 0.d0
        Esz = 0.d0
        Eex = 0.d0
        Eey = 0.d0
        Eez = 0.d0
        Bx = Bsx + Bex
        By = Bsy + Bey
        Bz = Bsz + Bez
        Ex = Esx + Eex
        Ey = Esy + Eey
        Ez = Esz + Eez
        Te = 1.d0
        Ti = 1.d0
        eta = 0.d0
        call current
        call pressure

        ! controling parameter
        is_abnormal_resistance = .false.
        
    end subroutine



    subroutine stepon
        implicit none

        call continuity_equation
        call momentum_equation
        call energy_equation
        call EMField_eqution

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
        K2 = -tau * (vex * central_difference_x(Te+K1/2.0d0,hx) + vey * central_difference_y(Te+K1/2.0d0,hy) + vez * central_difference_z(Te+K1/2.0d0,hz))
        K3 = -tau * (vex * central_difference_x(Te+K2/2.0d0,hx) + vey * central_difference_y(Te+K2/2.0d0,hy) + vez * central_difference_z(Te+K2/2.0d0,hz))
        K4 = -tau * (vex * central_difference_x(Te+K3,hx) + vey * central_difference_y(Te+K3,hy) + vez * central_difference_z(Te+K3,hz))
        temp = Te + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(Te,hx)
        cdiffy = central_difference_y(Te,hy)
        cdiffz = central_difference_z(Te,hz)

        temp(2:4,:,:) = Te(2:4,:,:) - tau * (vex(2:4,:,:) * cdiffx(2:4,:,:) + vey(2:4,:,:) * cdiffy(2:4,:,:) + vez(2:4,:,:) * cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = Te(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:) * cdiffx(Nx-3:Nx-1,:,:) + vey(Nx-3:Nx-1,:,:) * cdiffy(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:) * cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = Te(:,2:4,:) - tau * (vex(:,2:4,:) * cdiffx(:,2:4,:) + vey(:,2:4,:) * cdiffy(:,2:4,:) + vez(:,2:4,:) * cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = Te(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:) * cdiffx(:,Ny-3:Ny-1,:) + vey(:,Ny-3:Ny-1,:) * cdiffy(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:) * cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = Te(:,:,2:4) - tau * (vex(:,:,2:4) * cdiffx(:,:,2:4) + vey(:,:,2:4) * cdiffy(:,:,2:4) + vez(:,:,2:4) * cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = Te(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1) * cdiffx(:,:,Nz-3:Nz-1) + vey(:,:,Nz-3:Nz-1) * cdiffy(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1) * cdiffz(:,:,Nz-3:Nz-1))

        ! add other terms
        temp = temp + tau * 2.0d0 / (3.0d0*ne) * (-pre * (central_difference_x(vex,hx) + central_difference_y(vey,hy) + central_difference_z(vez,hz)))

        ! boundary points
        call boundary(temp)

        ! update Te
        Te = temp

        ! Ti second
        ! inner points use RK4
        K1 = -tau * (vix * central_difference_x(Ti,hx) + viy * central_difference_y(Ti,hy) + viz * central_difference_z(Ti,hz))
        K2 = -tau * (vix * central_difference_x(Ti+K1/2.0d0,hx) + viy * central_difference_y(Ti+K1/2.0d0,hy) + viz * central_difference_z(Ti+K1/2.0d0,hz))
        K3 = -tau * (vix * central_difference_x(Ti+K2/2.0d0,hx) + viy * central_difference_y(Ti+K2/2.0d0,hy) + viz * central_difference_z(Ti+K2/2.0d0,hz))
        K4 = -tau * (vix * central_difference_x(Ti+K3,hx) + viy * central_difference_y(Ti+K3,hy) + vez * central_difference_z(Te+K3,hz))
        temp = Ti + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! subouter points use Euler
        cdiffx = central_difference_x(Ti,hx)
        cdiffy = central_difference_y(Ti,hy)
        cdiffz = central_difference_z(Ti,hz)

        temp(2:4,:,:) = Ti(2:4,:,:) - tau * (vix(2:4,:,:) * cdiffx(2:4,:,:) + viy(2:4,:,:) * cdiffy(2:4,:,:) + viz(2:4,:,:) * cdiffz(2:4,:,:))
        temp(Nx-3:Nx-1,:,:) = Ti(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:) * cdiffx(Nx-3:Nx-1,:,:) + viy(Nx-3:Nx-1,:,:) * cdiffy(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:) * cdiffz(Nx-3:Nx-1,:,:))
        temp(:,2:4,:) = Ti(:,2:4,:) - tau * (vix(:,2:4,:) * cdiffx(:,2:4,:) + viy(:,2:4,:) * cdiffy(:,2:4,:) + viz(:,2:4,:) * cdiffz(:,2:4,:))
        temp(:,Ny-3:Ny-1,:) = Ti(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:) * cdiffx(:,Ny-3:Ny-1,:) + viy(:,Ny-3:Ny-1,:) * cdiffy(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:) * cdiffz(:,Ny-3:Ny-1,:))
        temp(:,:,2:4) = Ti(:,:,2:4) - tau * (vix(:,:,2:4) * cdiffx(:,:,2:4) + viy(:,:,2:4) * cdiffy(:,:,2:4) + viz(:,:,2:4) * cdiffz(:,:,2:4))
        temp(:,:,Nz-3:Nz-1) = Te(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1) * cdiffx(:,:,Nz-3:Nz-1) + viy(:,:,Nz-3:Nz-1) * cdiffy(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1) * cdiffz(:,:,Nz-3:Nz-1))

        ! add other terms
        temp = temp + tau * 2.0d0 / (3.0d0*ni) * (-pri * (central_difference_x(vix,hx) + central_difference_y(viy,hy) + central_difference_z(viz,hz)))

        ! boundary points
        call boundary(temp)

        ! update Te
        Te = temp

        deallocate(K1, K2, K3, K4)
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
        K2 = -tau * (central_difference_x((ne+K1/2.0d0)*vex, hx)+central_difference_y((ne+K1/2.0d0)*vey, hy)+central_difference_z((ne+K1/2.0d0)*vez, hz))
        K3 = -tau * (central_difference_x((ne+K2/2.0d0)*vex, hx)+central_difference_y((ne+K2/2.0d0)*vey, hy)+central_difference_z((ne+K2/2.0d0)*vez, hz))
        K4 = -tau * (central_difference_x((ne+K3)*vex, hx)+central_difference_y((ne+K3)*vey, hy)+central_difference_z((ne+K3)*vez, hz))
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
        K2 = -tau * (central_difference_x((ni+K1/2.0d0)*vix, hx)+central_difference_y((ni+K1/2.0d0)*viy, hy)+central_difference_z((ni+K1/2.0d0)*viz, hz))
        K3 = -tau * (central_difference_x((ni+K2/2.0d0)*vix, hx)+central_difference_y((ni+K2/2.0d0)*viy, hy)+central_difference_z((ni+K2/2.0d0)*viz, hz))
        K4 = -tau * (central_difference_x((ni+K3)*vix, hx)+central_difference_y((ni+K3)*viy, hy)+central_difference_z((ni+K3)*viz, hz))
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
        temp(:,:,Nz-3:Nz-1) = ne(:,:,Nz-3:Nz-1) - tau * (cdiffx(:,:,Nz-3:Nz-1) + cdiffy(:,:,Nz-3:Nz-1) + cdiffz(:,:,Nz-3:Nz-1))
        ! boundary points
        call boundary(temp)

        ! update ni
        ni = temp

        deallocate(K1, K2, K3, K4, temp)
    end subroutine

    subroutine momentum_equation
        implicit none
        real(kind=8), allocatable :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
        real(kind=8), allocatable :: K1(:,:,:), K2(:,:,:), K3(:,:,:), K4(:,:,:)

        allocate(temp1(Nx,Ny,Nz), temp2(Nx,Ny,Nz), temp3(Nx,Ny,Nz))
        allocate(K1(Nx,Ny,Nz), K2(Nx,Ny,Nz), K3(Nx,Ny,Nz), K4(Nx,Ny,Nz))

        ! electron
        ! inner points use RK4
        ! vex
        K1 = -tau * (vex * central_difference_x(vex, hx) + vey * central_difference_y(vex, hy) + vez * central_difference_z(vex, hz))
        K2 = -tau * ((vex+K1/2.0d0) * central_difference_x(vex+K1/2.0d0, hx) + vey * central_difference_y(vex+K1/2.0d0, hy) + vez * central_difference_z(vex+K1/2.0d0, hz))
        K3 = -tau * ((vex+K2/2.0d0) * central_difference_x(vex+K2/2.0d0, hx) + vey * central_difference_y(vex+K2/2.0d0, hy) + vez * central_difference_z(vex+K2/2.0d0, hz))
        K4 = -tau * ((vex+K3) * central_difference_x(vex+K3, hx) + vey * central_difference_y(vex+K3, hy) + vez * central_difference_z(vex+K3, hz))
        temp1 = vex + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! vey
        K1 = -tau * (vex * central_difference_x(vey, hx) + vey * central_difference_y(vey, hy) + vez * central_difference_z(vey, hz))
        K2 = -tau * (vex * central_difference_x(vey+K1/2.0d0, hx) + (vey+K1/2.0d0) * central_difference_y(vey+K1/2.0d0, hy) + vez * central_difference_z(vey+K1/2.0d0, hz))
        K3 = -tau * (vex * central_difference_x(vey+K2/2.0d0, hx) + (vey+K2/2.0d0) * central_difference_y(vey+K2/2.0d0, hy) + vez * central_difference_z(vey+K2/2.0d0, hz))
        K4 = -tau * (vex * central_difference_x(vey+K3, hx) + (vey+K3) * central_difference_y(vey+K3, hy) + vez * central_difference_z(vey+K3, hz))
        temp2 = vey + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! vez
        K1 = -tau * (vex * central_difference_x(vez, hx) + vey * central_difference_y(vez, hy) + vez * central_difference_z(vez, hz))
        K2 = -tau * (vex * central_difference_x(vez+K1/2.0d0, hx) + vey * central_difference_y(vez+K1/2.0d0, hy) + (vez+K1/2.0d0) * central_difference_z(vez+K1/2.0d0, hz))
        K3 = -tau * (vex * central_difference_x(vez+K2/2.0d0, hx) + vey * central_difference_y(vez+K2/2.0d0, hy) + (vez+K2/2.0d0) * central_difference_z(vez+K2/2.0d0, hz))
        K4 = -tau * (vex * central_difference_x(vez+K3, hx) + vey * central_difference_y(vez+K3, hy) + (vez+K3) * central_difference_z(vez+K3, hz))
        temp3 = vez + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        

        ! subouter points use Euler
        call momentum_subouter(temp1, temp2, temp3, 1)

        ! add other terms
        temp1 = temp1 + tau / (me * ne) * (ne * qe * (Ex + (vey * Bz - vez * By)) - central_difference_x(pre, hx))
        temp2 = temp2 + tau / (me * ne) * (ne * qe * (Ey + (vez * Bx - vex * Bz)) - central_difference_y(pre, hy))
        temp3 = temp2 + tau / (me * ne) * (ne * qe * (Ez + (vex * By - vey * Bx)) - central_difference_z(pre, hz))

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
        K1 = -tau * (vix * central_difference_x(vix, hx) + viy * central_difference_y(vix, hy) + viz * central_difference_z(vix, hz))
        K2 = -tau * ((vix+K1/2.0d0) * central_difference_x(vix+K1/2.0d0, hx) + viy * central_difference_y(vix+K1/2.0d0, hy) + viz * central_difference_z(vix+K1/2.0d0, hz))
        K3 = -tau * ((vix+K2/2.0d0) * central_difference_x(vix+K2/2.0d0, hx) + viy * central_difference_y(vix+K2/2.0d0, hy) + viz * central_difference_z(vix+K2/2.0d0, hz))
        K4 = -tau * ((vix+K3) * central_difference_x(vix+K3, hx) + viy * central_difference_y(vix+K3, hy) + viz * central_difference_z(vix+K3, hz))
        temp1 = vix + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0

        ! viy
        K1 = -tau * (vix * central_difference_x(viy, hx) + viy * central_difference_y(viy, hy) + viz * central_difference_z(viy, hz))
        K2 = -tau * (vix * central_difference_x(viy+K1/2.0d0, hx) + (viy+K1/2.0d0) * central_difference_y(viy+K1/2.0d0, hy) + viz * central_difference_z(viy+K1/2.0d0, hz))
        K3 = -tau * (vix * central_difference_x(viy+K2/2.0d0, hx) + (viy+K2/2.0d0) * central_difference_y(viy+K2/2.0d0, hy) + viz * central_difference_z(viy+K2/2.0d0, hz))
        K4 = -tau * (vix * central_difference_x(viy+K3, hx) + (viy+K3) * central_difference_y(viy+K3, hy) + viz * central_difference_z(viy+K3, hz))
        temp2 = viy + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        
        ! viz
        K1 = -tau * (vix * central_difference_x(viz, hx) + viy * central_difference_y(viz, hy) + viz * central_difference_z(viz, hz))
        K2 = -tau * (vix * central_difference_x(viz+K1/2.0d0, hx) + viy * central_difference_y(viz+K1/2.0d0, hy) + (viz+K1/2.0d0) * central_difference_z(viz+K1/2.0d0, hz))
        K3 = -tau * (vix * central_difference_x(viz+K2/2.0d0, hx) + viy * central_difference_y(viz+K2/2.0d0, hy) + (viz+K2/2.0d0) * central_difference_z(viz+K2/2.0d0, hz))
        K4 = -tau * (vix * central_difference_x(viz+K3, hx) + viy * central_difference_y(viz+K3, hy) + (viz+K3) * central_difference_z(viz+K3, hz))
        temp3 = viz + (K1 + 2.0d0*K2 + 2.0d0*K3 + K4) / 6.0d0
        

        ! subouter points use Euler
        call momentum_subouter(temp1, temp2, temp3, 2)

        ! add other terms
        temp1 = temp1 + tau / (mi * ni) * (ni * qi * (Ex + (viy * Bz - viz * By)) - central_difference_x(pri, hx))
        temp2 = temp2 + tau / (mi * ni) * (ni * qi * (Ey + (viz * Bx - vix * Bz)) - central_difference_y(pri, hy))
        temp3 = temp3 + tau / (mi * ni) * (ni * qi * (Ez + (vix * By - viy * Bx)) - central_difference_z(pri, hz))

        ! boundary points
        call boundary(temp1)
        call boundary(temp2)
        call boundary(temp3)

        ! update
        vix = temp1
        viy = temp2
        viz = temp3

        deallocate(temp1, temp2, temp3)
        deallocate(K1, K2, K3, K4)
    end subroutine

    subroutine momentum_subouter(vx, vy, vz, species)
        implicit none
        real(kind=8) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        integer :: species  ! 1 for electron, 2 for ion
        real(kind=8), allocatable :: cdiffxx(:,:,:), cdiffxy(:,:,:), cdiffxz(:,:,:)
        real(kind=8), allocatable :: cdiffyx(:,:,:), cdiffyy(:,:,:), cdiffyz(:,:,:)
        real(kind=8), allocatable :: cdiffzx(:,:,:), cdiffzy(:,:,:), cdiffzz(:,:,:)
        integer :: size1, size2, size3

        allocate(cdiffxx(Nx,Ny,Nz), cdiffxy(Nx,Ny,Nz), cdiffxz(Nx,Ny,Nz))
        allocate(cdiffyx(Nx,Ny,Nz), cdiffyy(Nx,Ny,Nz), cdiffyz(Nx,Ny,Nz))
        allocate(cdiffzx(Nx,Ny,Nz), cdiffzy(Nx,Ny,Nz), cdiffzz(Nx,Ny,Nz))
        size1 = size(vx, 1)
        size2 = size(vx, 2)
        size3 = size(vx, 3)

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

            vx(2:4,:,:) = vex(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxx(2:4,:,:)+ vey(2:4,:,:)*cdiffyx(2:4,:,:) + vez(2:4,:,:)*cdiffzx(2:4,:,:))
            vx(Nx-3:Nx-1,:,:) = vex(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxx(Nx-3:Nx-1,:,:)+ vey(Nx-3:Nx-1,:,:)*cdiffyx(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzx(Nx-3:Nx-1,:,:))
            vx(:,2:4,:) = vex(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxx(:,2:4,:)+ vey(:,2:4,:)*cdiffyx(:,2:4,:) + vez(:,2:4,:)*cdiffzx(:,2:4,:))
            vx(:,Ny-3:Ny-1,:) = vex(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxx(:,Ny-3:Ny-1,:)+ vey(:,Ny-3:Ny-1,:)*cdiffyx(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzx(:,Ny-3:Ny-1,:))
            vx(:,:,2:4) = vex(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxx(:,:,2:4)+ vey(:,:,2:4)*cdiffyx(:,:,2:4) + vez(:,:,2:4)*cdiffzx(:,:,2:4))
            vx(:,:,Nz-3:Nz-1) = vex(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxx(:,:,Nz-3:Nz-1)+ vey(:,:,Nz-3:Nz-1)*cdiffyx(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzx(:,:,Nz-3:Nz-1))
            
            vy(2:4,:,:) = vey(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxy(2:4,:,:)+ vey(2:4,:,:)*cdiffyy(2:4,:,:) + vez(2:4,:,:)*cdiffzy(2:4,:,:))
            vy(Nx-3:Nx-1,:,:) = vey(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxy(Nx-3:Nx-1,:,:)+ vey(Nx-3:Nx-1,:,:)*cdiffyy(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzy(Nx-3:Nx-1,:,:))
            vy(:,2:4,:) = vey(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxy(:,2:4,:)+ vey(:,2:4,:)*cdiffyy(:,2:4,:) + vez(:,2:4,:)*cdiffzy(:,2:4,:))
            vy(:,Ny-3:Ny-1,:) = vey(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxy(:,Ny-3:Ny-1,:)+ vey(:,Ny-3:Ny-1,:)*cdiffyy(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzy(:,Ny-3:Ny-1,:))
            vy(:,:,2:4) = vey(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxy(:,:,2:4)+ vey(:,:,2:4)*cdiffyy(:,:,2:4) + vez(:,:,2:4)*cdiffzy(:,:,2:4))
            vy(:,:,Nz-3:Nz-1) = vey(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxy(:,:,Nz-3:Nz-1)+ vey(:,:,Nz-3:Nz-1)*cdiffyy(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzy(:,:,Nz-3:Nz-1))
            
            vz(2:4,:,:) = vez(2:4,:,:) - tau * (vex(2:4,:,:)*cdiffxz(2:4,:,:)+ vey(2:4,:,:)*cdiffyz(2:4,:,:) + vez(2:4,:,:)*cdiffzz(2:4,:,:))
            vz(Nx-3:Nx-1,:,:) = vez(Nx-3:Nx-1,:,:) - tau * (vex(Nx-3:Nx-1,:,:)*cdiffxz(Nx-3:Nx-1,:,:)+ vey(Nx-3:Nx-1,:,:)*cdiffyz(Nx-3:Nx-1,:,:) + vez(Nx-3:Nx-1,:,:)*cdiffzz(Nx-3:Nx-1,:,:))
            vz(:,2:4,:) = vez(:,2:4,:) - tau * (vex(:,2:4,:)*cdiffxz(:,2:4,:)+ vey(:,2:4,:)*cdiffyz(:,2:4,:) + vez(:,2:4,:)*cdiffzz(:,2:4,:))
            vz(:,Ny-3:Ny-1,:) = vez(:,Ny-3:Ny-1,:) - tau * (vex(:,Ny-3:Ny-1,:)*cdiffxz(:,Ny-3:Ny-1,:)+ vey(:,Ny-3:Ny-1,:)*cdiffyz(:,Ny-3:Ny-1,:) + vez(:,Ny-3:Ny-1,:)*cdiffzz(:,Ny-3:Ny-1,:))
            vz(:,:,2:4) = vez(:,:,2:4) - tau * (vex(:,:,2:4)*cdiffxz(:,:,2:4)+ vey(:,:,2:4)*cdiffyz(:,:,2:4) + vez(:,:,2:4)*cdiffzz(:,:,2:4))
            vz(:,:,Nz-3:Nz-1) = vez(:,:,Nz-3:Nz-1) - tau * (vex(:,:,Nz-3:Nz-1)*cdiffxz(:,:,Nz-3:Nz-1)+ vey(:,:,Nz-3:Nz-1)*cdiffyz(:,:,Nz-3:Nz-1) + vez(:,:,Nz-3:Nz-1)*cdiffzz(:,:,Nz-3:Nz-1))
   
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

            vx(2:4,:,:) = vix(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxx(2:4,:,:)+ viy(2:4,:,:)*cdiffyx(2:4,:,:) + viz(2:4,:,:)*cdiffzx(2:4,:,:))
            vx(Nx-3:Nx-1,:,:) = vix(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxx(Nx-3:Nx-1,:,:)+ viy(Nx-3:Nx-1,:,:)*cdiffyx(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzx(Nx-3:Nx-1,:,:))
            vx(:,2:4,:) = vix(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxx(:,2:4,:)+ viy(:,2:4,:)*cdiffyx(:,2:4,:) + viz(:,2:4,:)*cdiffzx(:,2:4,:))
            vx(:,Ny-3:Ny-1,:) = vix(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxx(:,Ny-3:Ny-1,:)+ viy(:,Ny-3:Ny-1,:)*cdiffyx(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzx(:,Ny-3:Ny-1,:))
            vx(:,:,2:4) = vix(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxx(:,:,2:4)+ viy(:,:,2:4)*cdiffyx(:,:,2:4) + viz(:,:,2:4)*cdiffzx(:,:,2:4))
            vx(:,:,Nz-3:Nz-1) = vix(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxx(:,:,Nz-3:Nz-1)+ viy(:,:,Nz-3:Nz-1)*cdiffyx(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzx(:,:,Nz-3:Nz-1))
            
            vy(2:4,:,:) = viy(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxy(2:4,:,:)+ viy(2:4,:,:)*cdiffyy(2:4,:,:) + viz(2:4,:,:)*cdiffzy(2:4,:,:))
            vy(Nx-3:Nx-1,:,:) = viy(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxy(Nx-3:Nx-1,:,:)+ viy(Nx-3:Nx-1,:,:)*cdiffyy(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzy(Nx-3:Nx-1,:,:))
            vy(:,2:4,:) = viy(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxy(:,2:4,:)+ viy(:,2:4,:)*cdiffyy(:,2:4,:) + viz(:,2:4,:)*cdiffzy(:,2:4,:))
            vy(:,Ny-3:Ny-1,:) = viy(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxy(:,Ny-3:Ny-1,:)+ viy(:,Ny-3:Ny-1,:)*cdiffyy(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzy(:,Ny-3:Ny-1,:))
            vy(:,:,2:4) = viy(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxy(:,:,2:4)+ viy(:,:,2:4)*cdiffyy(:,:,2:4) + viz(:,:,2:4)*cdiffzy(:,:,2:4))
            vy(:,:,Nz-3:Nz-1) = viy(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxy(:,:,Nz-3:Nz-1)+ viy(:,:,Nz-3:Nz-1)*cdiffyy(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzy(:,:,Nz-3:Nz-1))
            
            vz(2:4,:,:) = viz(2:4,:,:) - tau * (vix(2:4,:,:)*cdiffxz(2:4,:,:)+ viy(2:4,:,:)*cdiffyz(2:4,:,:) + viz(2:4,:,:)*cdiffzz(2:4,:,:))
            vz(Nx-3:Nx-1,:,:) = viz(Nx-3:Nx-1,:,:) - tau * (vix(Nx-3:Nx-1,:,:)*cdiffxz(Nx-3:Nx-1,:,:)+ viy(Nx-3:Nx-1,:,:)*cdiffyz(Nx-3:Nx-1,:,:) + viz(Nx-3:Nx-1,:,:)*cdiffzz(Nx-3:Nx-1,:,:))
            vz(:,2:4,:) = viz(:,2:4,:) - tau * (vix(:,2:4,:)*cdiffxz(:,2:4,:)+ viy(:,2:4,:)*cdiffyz(:,2:4,:) + viz(:,2:4,:)*cdiffzz(:,2:4,:))
            vz(:,Ny-3:Ny-1,:) = viz(:,Ny-3:Ny-1,:) - tau * (vix(:,Ny-3:Ny-1,:)*cdiffxz(:,Ny-3:Ny-1,:)+ viy(:,Ny-3:Ny-1,:)*cdiffyz(:,Ny-3:Ny-1,:) + viz(:,Ny-3:Ny-1,:)*cdiffzz(:,Ny-3:Ny-1,:))
            vz(:,:,2:4) = viz(:,:,2:4) - tau * (vix(:,:,2:4)*cdiffxz(:,:,2:4)+ viy(:,:,2:4)*cdiffyz(:,:,2:4) + viz(:,:,2:4)*cdiffzz(:,:,2:4))
            vz(:,:,Nz-3:Nz-1) = viz(:,:,Nz-3:Nz-1) - tau * (vix(:,:,Nz-3:Nz-1)*cdiffxz(:,:,Nz-3:Nz-1)+ viy(:,:,Nz-3:Nz-1)*cdiffyz(:,:,Nz-3:Nz-1) + viz(:,:,Nz-3:Nz-1)*cdiffzz(:,:,Nz-3:Nz-1))

        end if
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

        print *, divB_max
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

        real(kind=8), allocatable ::  temp1(:,:,:), temp2(:,:,:)
        real(kind=8) :: dxyz, dt1, dt2

        allocate(temp1(Nx,Ny,Nz), temp2(Nx,Ny,Nz))
        !call foreta(time, 1)
        call pressure

        dxyz = .5*hx*hy*hz/sqrt((hx**2*hy**2 + hy**2*hz**2 + hx**2*hz**2))

        temp1 = dxyz / (sqrt(vex**2 + vey**2 + vez**2) + sqrt(Bx**2 + By**2 + Bz**2 + gamma*pre)/(me*ne))
        temp2 = dxyz / (sqrt(vix**2 + viy**2 + viz**2) + sqrt(Bx**2 + By**2 + Bz**2 + gamma*pri)/(mi*ni))

        dt1 = minval(minval(minval(temp1,3),2),1)
        dt2 = minval(minval(minval(temp2,3),2),1)

        tau = 0.5*min(dt1, dt2)
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

    subroutine current
        ! calculate the current
        implicit none

        jx = const1 * (ne*qe*vex + ni*qi*vix)
        jy = const1 * (ne*qe*vey + ni*qi*viy)
        jz = const1 * (ne*qe*vez + ni*qi*viz)
        
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

    subroutine grid

        implicit none

        integer :: jx, jy, jz

        hx = (xmax - xmin)/(Nx-1)
        hy = (ymax - ymin)/(Ny-1)
        hz = (zmax - zmin)/(Nz-1)

        do jx = 1, Nx
            x(jx) = xmin + (jx - 1)*hx
        end do

        do jy = 1, Ny
            y(jy) = ymin + (jy - 1)*hy
        end do

        do jz = 1, Nz
            z(jz) = zmin + (jz - 1)*hz
        end do

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

    subroutine pressure
        ! calculate the pressure of electron and ion
        implicit none

        pre = ne * Te
        pri = ni * Ti
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

    subroutine abnormal_resistance
        implicit none

        if(.not. is_abnormal_resistance) return

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
