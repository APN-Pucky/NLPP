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
        
    subroutine io_selectLoop(array)
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

    function io_getFileName(default_fname) result(fname)
        character(len=*), optional :: default_fname
	    character(len=100) :: fname
	    if(.NOT. present(default_fname))default_fname="" 
        write(*,'(3(A))', advance='no') "Dateiname (default=", default_fname ,"):"       
        read(*,'(100A)') fname
        if(fname=='') fname = default_fname
        write(*,'(2(A))') "Öffne: ", fname
    end function
    
end module
