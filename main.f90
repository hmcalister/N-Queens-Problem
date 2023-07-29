program NQueenProblem
    use iso_fortran_env, only: real64
    use calculations
    implicit none

    integer, parameter :: N = 100
    character(len=5) :: StringN
    character(len=20) :: matrixPrintFormatString
    real(real64), dimension(N, N) :: V
    type(NQueenParameters) :: params
    integer :: i
    integer :: j
    real(real64) :: sum

    params%dim=N
    params%A=0.5
    params%B=0.5
    params%C=0.5
    params%N_Neg=real(N, kind=real64)

    write(StringN, '(I4)') N
    matrixPrintFormatString = "("//StringN//"(F7.4, 2X))"

    V(:, :) = 0
    
    ! print matrixPrintFormatString, V
    
    sum = 0
    do i = 1, params%dim
        do j = 1, params%dim
            sum = sum + params%unit_activation(V, i, j)
        end do
    end do

    print *, sum/N**2

end program NQueenProblem