program main
    use io
    use leq
    implicit none
    call io_selectLoop((/&
        Select("TESTE LEQ", test_leq), &
        Select("NULL", test_leq) &
    /))

    contains 

    subroutine b()
        write(*,*) "b"
    end subroutine
    subroutine a()
        write(*,*) "a"
    end subroutine
end program
