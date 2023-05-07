program main
    !             3-D COMPRESSIBLE MHD CODE
    !    ************************************************
    !
    ! 1. Scheme  ---Modify 2 step Lax-Wendroff.
    ! 2. History --- completed in June, 1991.
    ! 3. Author  --- Zhiwei Ma (Geophysical Institute, UAF, AK 99775)
    ! 4. Project ---- Asymmetric Local magnetic reconnection.
    !```````````````````````````````````````````````````````````
    ! 5. revised by cjxiao on 2004-04-20 (in subroutine init):
    !    to chang X-line not perpendicular to MR plane.
    !
    !
    !                  Introduction
    ! The present version of this code simulate one half of the
    ! simulation box in x-y-z direction. In the six boundary,two boundaries
    ! at x=Lx and x=-Lx use inflow ones,three at y=Ly, y=-Ly , and z=Lz
    ! free (outgoing) ones,and one at z=-Lz symmetric (or antisymmetric )
    ! ones.
    !
    !     **********************************************
    !-----------------
    !     Basic variable definitions:
    !     x(mx,my,mz,1) to x(mx,my,mz,8) represents, respectively,
    !     rho, rho*vx, rho*vy, rho*vz, bx, by, bz and energy.
    !-----------------
    
    ! usage : gfortran [-On] Define_variables.f90 Custom_functions.f90 Custom_subroutines.f90 Two_Fluid_MHD.f90 -o MHD
    
    use define_variables
    use custom_functions
    use custom_subroutines

    implicit none

    integer :: runtime_start, runtime_end

    call system_clock(runtime_start)

    call initialize

    do while (nstep < nmax)
        call setdt
        call stepon

        nstep = nstep + 1
        time = time + tau
        print *, 'nstep= ', nstep, ' ', 'time= ', time, '', 'dt=', tau
    end do

    call destroy

    call system_clock(runtime_end)
    print *, 'runtime= ', (runtime_end-runtime_start)/1000.0, 's'

end program