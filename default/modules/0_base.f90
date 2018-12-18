module base
    implicit none

    interface assert
        module procedure assert_int_ab, assert_int_abs
    end interface
 
    contains

    subroutine assert_int_abs(a,b,s)
        integer :: a,b
        character(len=*) :: s
        if(a/=b) then
            print *, s , a , "/=" , b
            call EXIT(1)
        end if
    end subroutine

    subroutine assert_int_ab(a,b)
        integer :: a,b
        call assert(a,b,"Assert failed: ")
    end subroutine

    subroutine error(s)
        character(len=*) :: s
        print *,s
        call exit(1)
    end subroutine


end module
