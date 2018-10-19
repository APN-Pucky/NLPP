program main
    use io
    implicit none
    call select_loop((/&
        Select("Test", a),&
        Select("BBB", b)&
    /))

    contains 

    subroutine b()
        write(*,*) "b"
    end subroutine
    subroutine a()
        write(*,*) "a"
    end subroutine
end program
