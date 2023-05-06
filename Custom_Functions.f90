module custom_functions
    !implicit none

contains

    character*3 function cn(n)
        !-----assume that n is no greater than 999
        !
        !
        !-----separate the digits
        implicit none
        integer :: n, n1, n2, n3

        n1 = n/100
        n2 = (n - 100*n1)/10
        n3 = n - 100*n1 - 10*n2

        !-----stick together cn using char function

        n1 = n1 + 48
        n2 = n2 + 48
        n3 = n3 + 48
        cn(1:1) = char(n1)
        cn(2:2) = char(n2)
        cn(3:3) = char(n3)

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


end module Custom_functions
