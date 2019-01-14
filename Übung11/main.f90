program main
    use io
    use matrix
    use mfunktion
    implicit none
    call init
    call io_selectLoop((/&
        Select("Test a", test_a),&
        Select("Test b", test_b),&
        Select("Test DFT", test_a)&
    /))

    contains 

    subroutine test_a()
     
    end subroutine

    subroutine test_b()
        type(funktion) :: p
        integer :: i
        p = 1-heaviside(x**2-1)
        do i=1,10
            print *,real(i)/5, p%get(real(i)/5)
        end do
    end subroutine
end program
