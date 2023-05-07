program test
    implicit none
    integer :: i, j
    real(kind=8) :: x(2,2), y(2,2), z(2,2)

    x = reshape([1.0, 2.0, 3.0, 4.0], [2,2])
    y = reshape([2.0, 3.0, 4.0, 5.0], [2,2])
    z = reshape([3.0, 4.0, 5.0, 6.0], [2,2])

    write(*,"(3(1x,D9.2))") ((x(i,j), y(i,j), z(i,j), i=1,2), j=1,2)
end program test