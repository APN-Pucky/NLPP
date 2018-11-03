module IO 
    implicit none

    abstract interface
        subroutine sub_interface()
        end subroutine
    end interface
    type :: Select
        Character (len=20) :: name
        procedure(sub_interface),POINTER,NOPASS :: exec
    end type
    logical, dimension(99) :: nunit=.true.
    !formats!
    !100 FORMAT(50("="),T8,A)
    !101 FORMAT(50("-"),T5,A)
    contains
        
    subroutine io_selectLoop(array)
        type(Select), intent(in),dimension(:) :: array
        integer :: s
        s=-1
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

    function io_getInteger(default_number) result(num)
        integer, optional :: default_number
        integer :: num, st
        if(.NOT. present(default_number))default_number=0 
        write(*,*) "Integer (default=", default_number ,"):"       
        read(*,*, iostat=st) num
        if(st .GT. 0) num = default_number
        write(*,*) "Using: ", num
    end function
    
 
    function io_getFileName(default_fname) result(fname)
        character(len=*), optional :: default_fname
        character(len=100) :: fname
        if(.NOT. present(default_fname))default_fname="" 
        write(*,'(3(A))', advance='no') "File (default=", default_fname ,"):"       
        read(*,'(100A)') fname
        if(fname=='') fname = default_fname
        write(*,'(2(A))') "Opening: ", fname
    end function
    
    function io_openFile(default_fname,mode) result(i)
        integer ::i
        Character(len=*) :: mode,default_fname
        do i=1,99
            if(nunit(i)) then  
                nunit(i) = .false.
                open(unit=i, file=io_getFileName(default_fname), status=mode) 
                return 
            endif 
        enddo
        i= -1
    end function

    subroutine io_close(i)
        integer :: i
        nunit(i) = .true.
        close(i)
    end subroutine

    subroutine io_skip(f,n)
        integer :: n,j,f
        do j=1,n
            read(f,*) 
        end do
   end subroutine 
end module
