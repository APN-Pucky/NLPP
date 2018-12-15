program main
    use io
    implicit none
    call io_selectLoop((/&
        Select("Test a", a),&
        Select("Test b", b),&
        Select("Test c", c),&
    /))

    contains 

    subroutine b()
        call assert(1,2)
        write(*,*) "b"
    end subroutine
    subroutine a()
        integer :: ff
        real :: xyz1(3) = (/ 1,2,4/)
        real :: xyz2(3) = (/ 4,2,4/)
        write(*,*) "a"
        ff = io_openFile("data/test.xyz","replace")
        CALL io_writeXYZ(ff, (/ 'O','C'/), (/ xyz1,xyz2 /))
    end subroutine
end program
