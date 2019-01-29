module base
    implicit none

    interface assert
        module procedure assert_int_ab, assert_int_ab_msg, assert_logical_msg,assert_logical
    end interface
 
    contains

   function rtoa(r)
        real :: r
        character(len=1024) ::rtoa 
        write(rtoa,*) r
        rtoa = adjustl(rtoa)
    end function
        
    function itoa(i)
        integer :: i
        character(len=1024) ::itoa 
        write(itoa,*) i
        itoa = adjustl(itoa)
    end function

    subroutine assert_logical_msg(l,s)
        logical :: l
        character(len=*) :: s
        if(.NOT. l) then
            call error(s)
        end if
    end subroutine
 
    subroutine assert_int_ab_msg(a,b,s)
        integer :: a,b
        character(len=*) :: s
        call assert(a==b, s//": "//trim(itoa(a))//"/="//trim(itoa(b)))
    end subroutine

    subroutine assert_int_ab(a,b)
        integer :: a,b
        call assert(a,b,"Assert failed")
    end subroutine
    
    subroutine assert_logical(l)
        logical :: l
        call assert(l,"Assert failed")
    end subroutine

    subroutine error(s)
        character(len=*) :: s
        print *,s
        call exit(1)
    end subroutine


end module
