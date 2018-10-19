module io
    implicit none

    abstract interface
        subroutine sub_interface()
        end subroutine
    end interface
    type :: Select
        Character (len=20) :: name
        procedure(sub_interface),POINTER,NOPASS :: exec
    end type
    !formats!
    !100 FORMAT(50("="),T8,A)
    !101 FORMAT(50("-"),T5,A)
    contains
        
    subroutine select_loop(array)
        type(Select), intent(in),dimension(:) :: array
        integer :: s=-1
        do while(s /= 0)
            if(.NOT. (s> 0 .AND. s <= size(array))) then
                call print_opt(array%name)
            end if
            write(*,'(A)',advance="no") "> "
            read(*,*) s
            if(s> 0 .AND. s <= size(array)) then
                call array(s)%exec
            end if
        end do
        write(*,'(A)') "Quit"
    end subroutine

    subroutine print_opt(names)
        Character(len=20), dimension(:) :: names
        integer :: j
        Character(*),parameter :: f='(T2,I2,") ",A)'
        j=0
        write(*,'(A)') "Options: "
        write(*,f) j,"Quit"
        do j=1,size(names)
            write(*,f) j, names(j)
        end do
    end subroutine
end module
