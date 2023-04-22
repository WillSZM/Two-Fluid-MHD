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
    
    ! usage : gfortran [-On] Define_variables.f90 Custom_functions.f90 Custom_subroutines.f90 MHD.f90 -o MHD
    
    use define_variables
    use custom_functions
    use custom_subroutines

    implicit none

    call set_values

    ! call input

    call gridpnt
    
    call initia
    
    if (lrstrt) then
        call readin(nst2, cont, dtime)
    end if
 
    if (time .eq. 0) then
        call current(x,1)
        call recrd1
        !call facur
    end if

    do while (nst < nend)
        call setdt
        call stepon

        nstep = nstep + 1
        time = time + dt
        write (*, *) nstep, ' ', 'time=', time, '', 'dt=', dt         !zxg
        write (*, *) 'nst=', nst
        
        open (unit=16, file='stepnm', status='unknown', form='formatted')
        write (16, 99) nstep, time
99      format(i5, f9.5)
        close (16)
        !
        !zxg to continue
        if (abs(mod(time, 2.d0)) .le. dt) then
            open (unit=17, file='continue', status="unknown", form="unformatted")
            write (17) ncase, nstep, time, nst
            write (17) x
            close (17)
        end if
        !zxg continue end

        if (abs(time - nst*dtime) .le. dt .and. &
            (nstep - nstp(nst)) .ge. dnstep) then
            nst = nst + nint
            nstp(nst) = nstep
            call current(x,1)
            call recrd1
        end if

        !yg----------------
        if ((time .gt. t0) .and. ((time - t0) .le. dt)) then
            call incident_plasma(x, xi, time, t0, 0)
        end if

        if ((time .gt. t0) .and. ((time - t0) .gt. dt)) then
            call incident_plasma(x, xi, time, t0, 1)
        end if
        !yg---------------

        if (nstep .gt. nstop) exit
    end do

end program
