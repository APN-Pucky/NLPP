program main
    use io
    use verlet
    implicit none
    call io_selectLoop((/&
        Select("Test a", a),&
        Select("Test b", b),&
        Select("CO cms/lab", test_CO)&
    /))

    contains 
    subroutine b()
        call assert(1,1)
        write(*,*) "Test b"
    end subroutine
    subroutine a()
        integer :: ff
        real :: xyz1(3) = (/ 1,2,4/)
        real :: xyz2(3) = (/ 4,2,4/)
        write(*,*) "Test a"
        ff = io_openFile("data/test.xyz","replace")
        CALL io_writeXYZ(ff, (/ 'O','C'/), (/ xyz1,xyz2 /))
        call io_closeFile(ff)
    end subroutine
end program
