module custom_functions
    !implicit none

contains

    function int2str(num)
        ! assume that n is no greater than 999999
        implicit none

        integer :: num
        character(len=6) :: int2str

        write(int2str, "(I6.6)") num

    end function


    character*6 function cn(n)
        !-----assume that n is no greater than 999999
        !
        !
        !-----separate the digits
        implicit none
        integer :: n, n1, n2, n3, n4, n5, n6

        n1 = n/100000
        n2 = (n - 100000*n1)/10000
        n3 = (n - 100000*n1 - 10000*n2)/1000
        n4 = (n - 100000*n1 - 10000*n2 - 1000*n3)/100
        n5 = (n - 100000*n1 - 10000*n2 - 1000*n3 - 100*n4)/10
        n6 = n - 100000*n1 - 10000*n2 - 1000*n3 - 100*n4 - 10*n5

        !-----stick together cn using char function

        n1 = n1 + 48
        n2 = n2 + 48
        n3 = n3 + 48
        n4 = n4 + 48
        n5 = n5 + 48
        n6 = n6 + 48

        cn(1:1) = char(n1)
        cn(2:2) = char(n2)
        cn(3:3) = char(n3)
        cn(4:4) = char(n4)
        cn(5:5) = char(n5)
        cn(6:6) = char(n6)

        return
    end function

    function d1fc(fm, f0, fp, xm1, x0, xp1)
        implicit none
        real(kind=8) :: fm, f0, fp, xm1, x0, xp1
        real(kind=8) :: d1fc

        d1fc = ((xm1 - x0)/(xp1 - x0)*(fp - f0) - (xp1 - x0)/(xm1 - x0)*(fm - f0))/(xm1 - xp1)
    end function

    function central_difference_x(a,hx)
        !a is a 3D array, hx is the grid spacing in x direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: central_difference_x
        real(kind=8) :: hx

        central_difference_x(2:size(a,1)-1,:,:) = (a(3:size(a,1),:,:) - a(1:size(a,1)-2,:,:))/(2*hx)
        central_difference_x(1,:,:) = (a(2,:,:) - a(1,:,:))/hx
        central_difference_x(size(a,1),:,:) = (a(size(a,1),:,:) - a(size(a,1)-1,:,:))/hx
    end function

    function central_difference_y(a,hy)
        !a is a 3D array, hy is the grid spacing in y direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: central_difference_y
        real(kind=8) :: hy

        central_difference_y(:,2:size(a,2)-1,:) = (a(:,3:size(a,2),:) - a(:,1:size(a,2)-2,:))/(2*hy)
        central_difference_y(:,1,:) = (a(:,2,:) - a(:,1,:))/hy
        central_difference_y(:,size(a,2),:) = (a(:,size(a,2),:) - a(:,size(a,2)-1,:))/hy
    end function

    function central_difference_z(a,hz)
        !a is a 3D array, hz is the grid spacing in z direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: central_difference_z
        real(kind=8) :: hz

        central_difference_z(:,:,2:size(a,3)-1) = (a(:,:,3:size(a,3)) - a(:,:,1:size(a,3)-2))/(2*hz)
        central_difference_z(:,:,1) = (a(:,:,2) - a(:,:,1))/hz
        central_difference_z(:,:,size(a,3)) = (a(:,:,size(a,3)) - a(:,:,size(a,3)-1))/hz
    end function

    function cdiff2sd_x(a,hx)
        ! a is a 3D array, hx is the grid spacing in x direction
        ! central difference of second derivative in x direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff2sd_x
        real(kind=8) :: hx
        
        cdiff2sd_x(2:size(a,1)-1,:,:) = (a(3:size(a,1),:,:) - 2*a(2:size(a,1)-1,:,:) + a(1:size(a,1)-2,:,:))/(hx**2)
        cdiff2sd_x(1,:,:) = (-a(4,:,:) + 4*a(3,:,:) - 5*a(2,:,:) + 2*a(1,:,:))/(hx**2)
        cdiff2sd_x(size(a,1),:,:) = (a(size(a,1)-3,:,:) - 4*a(size(a,1)-2,:,:) + 5*a(size(a,1)-1,:,:) - 2*a(size(a,1),:,:))/(hx**2)
    end function

    function cdiff2sd_y(a,hy)
        ! a is a 3D array, hy is the grid spacing in y direction
        ! central difference of second derivative in y direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff2sd_y
        real(kind=8) :: hy
        
        cdiff2sd_y(:,2:size(a,2)-1,:) = (a(:,3:size(a,2),:) - 2*a(:,2:size(a,2)-1,:) + a(:,1:size(a,2)-2,:))/(hy**2)
        cdiff2sd_y(:,1,:) = (-a(:,4,:) + 4*a(:,3,:) - 5*a(:,2,:) + 2*a(:,1,:))/(hy**2)
        cdiff2sd_y(:,size(a,2),:) = (a(:,size(a,2)-3,:) - 4*a(:,size(a,2)-2,:) + 5*a(:,size(a,2)-1,:) - 2*a(:,size(a,2),:))/(hy**2)
    end function

    function cdiff2sd_z(a,hz)
        ! a is a 3D array, hz is the grid spacing in z direction
        ! central difference of second derivative in z direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff2sd_z
        real(kind=8) :: hz
        
        cdiff2sd_z(:,:,2:size(a,3)-1) = (a(:,:,3:size(a,3)) - 2*a(:,:,2:size(a,3)-1) + a(:,:,1:size(a,3)-2))/(hz**2)
        cdiff2sd_z(:,:,1) = (-a(:,:,4) + 4*a(:,:,3) - 5*a(:,:,2) + 2*a(:,:,1))/(hz**2)
        cdiff2sd_z(:,:,size(a,3)) = (a(:,:,size(a,3)-3) - 4*a(:,:,size(a,3)-2) + 5*a(:,:,size(a,3)-1) - 2*a(:,:,size(a,3)))/(hz**2)
    end function

    function cdiff4_x(a,hx)
        ! a is a 3D array, hx is the grid spacing in x direction
        ! 4th order central difference in x direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff4_x
        real(kind=8) :: hx

        cdiff4_x(3:size(a,1)-2,:,:) = (-a(5:size(a,1),:,:) + 8*a(4:size(a,1)-1,:,:) & 
                                        - 8*a(2:size(a,1)-3,:,:) + a(1:size(a,1)-4,:,:))/(12*hx)
        cdiff4_x(2,:,:) = (a(3,:,:) - a(1,:,:))/(2*hx)
        cdiff4_x(1,:,:) = (-3.0d0*a(1,:,:) + 4.0d0*a(2,:,:) - a(3,:,:))/(2*hx)
        cdiff4_x(size(a,1)-1,:,:) = (a(size(a,1),:,:) - a(size(a,1)-2,:,:))/(2*hx)
        cdiff4_x(size(a,1),:,:) = (3.0d0*a(size(a,1),:,:) - 4.0d0*a(size(a,1)-1,:,:) + a(size(a,1)-2,:,:))/(2*hx)
    end function

    function cdiff4_y(a,hy)
        ! a is a 3D array, hy is the grid spacing in y direction
        ! 4th order central difference in y direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff4_y
        real(kind=8) :: hy

        cdiff4_y(:,3:size(a,2)-2,:) = (-a(:,5:size(a,2),:) + 8*a(:,4:size(a,2)-1,:) & 
                                        - 8*a(:,2:size(a,2)-3,:) + a(:,1:size(a,2)-4,:))/(12*hy)
        cdiff4_y(:,2,:) = (a(:,3,:) - a(:,1,:))/(2*hy)
        cdiff4_y(:,1,:) = (-3.0d0*a(:,1,:) + 4.0d0*a(:,2,:) - a(:,3,:))/(2*hy)
        cdiff4_y(:,size(a,2)-1,:) = (a(:,size(a,2),:) - a(:,size(a,2)-2,:))/(2*hy)
        cdiff4_y(:,size(a,2),:) = (3.0d0*a(:,size(a,2),:) - 4.0d0*a(:,size(a,2)-1,:) + a(:,size(a,2)-2,:))/(2*hy)
    end function

    function cdiff4_z(a,hz)
        ! a is a 3D array, hz is the grid spacing in z direction
        ! 4th order central difference in z direction
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: a
        real(kind=8), dimension(size(a,1),size(a,2),size(a,3)) :: cdiff4_z
        real(kind=8) :: hz

        cdiff4_z(:,:,3:size(a,3)-2) = (-a(:,:,5:size(a,3)) + 8*a(:,:,4:size(a,3)-1) & 
                                        - 8*a(:,:,2:size(a,3)-3) + a(:,:,1:size(a,3)-4))/(12*hz)
        cdiff4_z(:,:,2) = (a(:,:,3) - a(:,:,1))/(2*hz)
        cdiff4_z(:,:,1) = (-3.0d0*a(:,:,1) + 4.0d0*a(:,:,2) - a(:,:,3))/(2*hz)
        cdiff4_z(:,:,size(a,3)-1) = (a(:,:,size(a,3)) - a(:,:,size(a,3)-2))/(2*hz)
        cdiff4_z(:,:,size(a,3)) = (3.0d0*a(:,:,size(a,3)) - 4.0d0*a(:,:,size(a,3)-1) + a(:,:,size(a,3)-2))/(2*hz)
    end function

    function Lax(x,f,g,h,s,tau,hx,hy,hz)
        ! x is unknown term, 3D array; f,g,h are flux terms; s is the source term
        ! Lax scheme
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: x,f,g,h,s
        real(kind=8), dimension(size(x,1),size(x,2),size(x,3)) :: Lax
        real(kind=8) :: tau,hx,hy,hz
        integer :: s1,s2,s3

        s1 = size(x,1)
        s2 = size(x,2)
        s3 = size(x,3)

        Lax(2:s1-1,2:s2-1,2:s3-1) = (x(3:s1,2:s2-1,2:s3-1) + x(1:s1-2,2:s2-1,2:s3-1) & 
                                    + x(2:s1-1,3:s2,2:s3-1) + x(2:s1-1,1:s2-2,2:s3-1) &
                                    + x(2:s1-1,2:s2-1,3:s3) + x(2:s1-1,2:s2-1,1:s3-2)) / 6.0d0 & 
                                    - tau/(2.0d0*hx)*(f(3:s1,2:s2-1,2:s3-1)-f(1:s1-2,2:s2-1,2:s3-1)) & 
                                    - tau/(2.0d0*hy)*(g(2:s1-1,3:s2,2:s3-1)-g(2:s1-1,1:s2-2,2:s3-1)) &
                                    - tau/(2.0d0*hz)*(h(2:s1-1,2:s2-1,3:s3)-h(2:s1-1,2:s2-1,1:s3-2)) &
                                    + tau*s(2:s1-1,2:s2-1,2:s3-1)
    end function

    function MacCormack(x,f,g,h,s,tau,hx,hy,hz)
        ! MacCormack scheme
        ! x is unknown term, 3D array; f,g,h are flux terms; s is the source term
        implicit none
        real(kind=8), dimension(:,:,:), intent(in) :: x,f,g,h,s
        real(kind=8), dimension(size(x,1),size(x,2),size(x,3)) :: MacCormack
        real(kind=8) :: tau,hx,hy,hz
        integer :: s1,s2,s3
        real(kind=8) :: a, b, c    ! parameter for MacCormack scheme
        real(kind=8), allocatable :: temp1(:,:,:), temp2(:,:,:)

        allocate(temp1(size(x,1),size(x,2),size(x,3)), temp2(size(x,1),size(x,2),size(x,3)))

        a = 1.0d0
        b = 1.0d0
        c = 1.0d0

        s1 = size(x,1)
        s2 = size(x,2)
        s3 = size(x,3)

        ! Prediction step
        temp1(2:s1-1,2:s2-1,2:s3-1) = x - tau/hx*(a*(f(3:s1,2:s2-1,2:s3-1)-f(2:s1-1,2:s2-1,2:s3-1)) & 
                                                    + (1-a)*(f(2:s1-1,2:s2-1,2:s3-1)-f(1:s1-2,2:s2-1,2:s3-1))) &
                                        - tau/hy*(b*(g(2:s1-1,3:s2,2:s3-1)-g(2:s1-1,2:s2-1,2:s3-1)) &
                                                    + (1-b)*(g(2:s1-1,2:s2-1,2:s3-1)-g(2:s1-1,1:s2-2,2:s3-1))) &
                                        - tau/hz*(c*(h(2:s1-1,2:s2-1,3:s3)-h(2:s1-1,2:s2-1,2:s3-1)) &
                                                    + (1-c)*(h(2:s1-1,2:s2-1,2:s3-1)-h(2:s1-1,2:s2-1,1:s3-2))) &
                                        + tau*s(2:s1-1,2:s2-1,2:s3-1)

        ! Correction step
        temp2(2:s1-1,2:s2-1,2:s3-1) = temp1(2:s1-1,2:s2-1,2:s3-1) &
                                            - tau/hx*((1-a)*(f(3:s1,2:s2-1,2:s3-1)-f(2:s1-1,2:s2-1,2:s3-1)) &
                                                        + a*(f(2:s1-1,2:s2-1,2:s3-1)-f(1:s1-2,2:s2-1,2:s3-1))) &
                                            - tau/hy*((1-b)*(g(2:s1-1,3:s2,2:s3-1)-g(2:s1-1,2:s2-1,2:s3-1)) &
                                                        + b*(g(2:s1-1,2:s2-1,2:s3-1)-g(2:s1-1,1:s2-2,2:s3-1))) &
                                            - tau/hz*((1-c)*(h(2:s1-1,2:s2-1,3:s3)-h(2:s1-1,2:s2-1,2:s3-1)) &
                                                        + c*(h(2:s1-1,2:s2-1,2:s3-1)-h(2:s1-1,2:s2-1,1:s3-2))) &
                                            + tau*s(2:s1-1,2:s2-1,2:s3-1)

        MacCormack(2:s1-1,2:s2-1,2:s3-1) = (temp1(2:s1-1,2:s2-1,2:s3-1) + temp2(2:s1-1,2:s2-1,2:s3-1)) / 2.0d0

        deallocate(temp1, temp2)
    end function
end module Custom_functions
