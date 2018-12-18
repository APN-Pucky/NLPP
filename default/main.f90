program main
    use io
    implicit none
    call init
    call io_selectLoop((/&
        Select("Test a", a),&
        Select("Test b", b)&
    /))

    contains 
    subroutine a()
        write(*,*) "a"
    end subroutine

    subroutine b()
        write(*,*) "b"
    end subroutine
end program
