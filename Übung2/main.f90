program main
    use io
    use metropolis
    implicit none
    call metro
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
