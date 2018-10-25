program main
    use mIO
    use metropolis
    use mMath
    implicit none
    call metro
    !call list_test
    !call math_test
    !call select_loop((/&
    !    Select("Elem Test",test),&
    !    Select("Test", a),&
    !    Select("BBB", b)&
    !/))

    contains 

    subroutine b()
        write(*,*) "b"
    end subroutine
    subroutine a()
        write(*,*) "a"
    end subroutine
end program
