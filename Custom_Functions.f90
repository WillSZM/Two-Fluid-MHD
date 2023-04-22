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

end module Custom_functions
