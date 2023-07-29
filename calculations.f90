module calculations
    use iso_fortran_env, only: real64
    implicit none
    
    type :: NQueenParameters
        integer :: dim
        real(real64) :: A
        real(real64) :: B
        real(real64) :: C
        real(real64) :: N_Neg

        contains
            procedure :: state_energy
            procedure :: unit_activation
    end type

contains

    function state_energy(self, matrix) result(energy)
        class(NQueenParameters), intent(in) :: self
        real(real64), dimension(:, :), intent(in) :: matrix
        real(real64) :: energy
        real(real64) :: currentTerm
        integer :: i
        integer :: j
        integer :: k

        energy = 0

        ! Term 1
        currentTerm = 0
        do i = 1 , self%dim
            do j = 1, self%dim
                do k = 1, self%dim
                    if (k==j) cycle
                    currentTerm = currentTerm + matrix(i,k)
                end do
                currentTerm = currentTerm * matrix(i,j)

                do k = 1, self%dim
                    if (k==i) cycle
                    currentTerm = currentTerm + matrix(k,j)
                end do
                currentTerm = currentTerm * matrix(i,j)
            end do
        end do
        currentTerm = currentTerm * self%A
        energy = energy + currentTerm

        ! TERM 2
        currentTerm = 0
        do i = 2 , self%dim
            do j = 1, i-1
                do k = i-j+1, self%dim
                    if (k==i) cycle
                    currentTerm = currentTerm + matrix(k, k-i+j)                    
                end do
                currentTerm = currentTerm * matrix(i,j)
            end do
        end do
        currentTerm = currentTerm * self%B
        energy = energy + currentTerm

        ! TERM 3
        currentTerm = 0
        do i = 1 , self%dim
            do j = i , self%dim
                do k = 1, self%dim+i-j
                    if (k==i) cycle
                    currentTerm = currentTerm + matrix(k, k-i+j)                    
                end do
                currentTerm = currentTerm * matrix(i,j)
            end do
        end do
        currentTerm = currentTerm * self%B
        energy = energy + currentTerm

        ! TERM 4
        currentTerm = 0
        do i = 1 , self%dim
            do j = self%dim-i+1 , self%dim
                do k = i+j-self%dim, self%dim
                    if (k==i) cycle
                    currentTerm = currentTerm + matrix(k, i+j-k)                    
                end do
                currentTerm = currentTerm * matrix(i,j)
            end do
        end do
        currentTerm = currentTerm * self%B
        energy = energy + currentTerm

        ! TERM 5
        currentTerm = 0
        do i = 1 , self%dim-1
            do j = 1, self%dim-i
                do k = 1, i+j-1
                    if (k==i) cycle
                    currentTerm = currentTerm + matrix(k, i+j-k)                    
                end do
                currentTerm = currentTerm * matrix(i,j)
            end do
        end do
        currentTerm = currentTerm * self%B
        energy = energy + currentTerm

        ! TERM 6
        currentTerm = 0
        do i = 1 , self%dim
            do j = 1, self%dim
                currentTerm = currentTerm + (matrix(i,j) - self%N_Neg)
            end do
        end do
        currentTerm = currentTerm**2 * self%C
        energy = energy + currentTerm

        energy = energy * 0.5
    end function state_energy

    function unit_activation(self, matrix, i, j) result(activation)
        class(NQueenParameters), intent(in) :: self
        real(real64), dimension(:, :), intent(in) :: matrix
        real(real64) :: activation
        integer :: i
        integer :: j

        real(real64) :: currentTerm
        integer :: k
        integer :: l

        currentTerm = 0
        do k = 1, self%dim 
            if (k==j) cycle
            currentTerm = currentTerm + matrix(i,k)
        end do
        do k = 1, self%dim 
            if (k==i) cycle
            currentTerm = currentTerm + matrix(k,j)
        end do
        currentTerm = currentTerm * self%A
        activation = activation + currentTerm

        currentTerm = 0
        if (i-j > 0) then
            do k = i-j+1, self%dim
                if (k==i) cycle
                currentTerm = currentTerm + matrix(k, k-i+j)
            end do
        else
            do k = 1, self%dim+i-j
                if (k==i) cycle
                currentTerm = currentTerm + matrix(k, k-i+j)
            end do
        end if
        currentTerm = currentTerm * self%B
        activation = activation + currentTerm

        currentTerm = 0
        if (i+j > self%dim) then
            do k = i+j-self%dim, self%dim
                if (k==i) cycle
                currentTerm = currentTerm + matrix(k, i+j-k)
            end do
        else
            do k = 1, i+j-1
                if (k==i) cycle
                currentTerm = currentTerm + matrix(k, i+j-k)
            end do
        end if
        currentTerm = currentTerm * self%B
        activation = activation + currentTerm

        currentTerm = 0
        do k = 1, self%dim
            do l = 1, self%dim
                currentTerm = currentTerm + (matrix(k,l) - self%N_Neg)
            end do
        end do
        currentTerm = currentTerm * self%C
        activation = activation + currentTerm

        activation = activation * (-1)            
    end function unit_activation

end module calculations