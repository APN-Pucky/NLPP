module io 
    use base
    implicit none

    abstract interface
        subroutine sub_interface()
        end subroutine
    end interface

    type :: Select
        Character (len=20) :: name
        procedure(sub_interface),POINTER,NOPASS :: exec
    end type

    integer, parameter :: max_n=64
    logical, dimension(max_n) :: nunit=.true.

    !formats!
    !100 FORMAT(50("="),T8,A)
    !101 FORMAT(50("-"),T5,A)
    contains
        
    subroutine io_selectLoop(array)
        type(Select), intent(in),dimension(:) :: array
        integer :: s,i
        character(*), parameter :: ff = '(50("-"),T8," ", A, " ")'
        character(*), parameter :: lf = '(50("="))'
        s=-1
        do while(s /= 99)
            if(.NOT. (s> 0 .AND. s <= size(array))) then
                call print_opt(array%name)
            end if
            write(*,'(A)',advance="no") "> "
            read(*,*) s

            if( s==0) then
                do i=1,size(array)
                    print lf
                    print ff, trim(array(i)%name)
                    print lf
                    call array(i)%exec
                end do
            end if

            if(s> 0 .AND. s <= size(array)) then
                print lf
                print ff, trim(array(s)%name)
                print lf
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
        write(*,f) 0,"All"
        write(*,f) 99,"Quit"
        do j=1,size(names)
            write(*,f) j, names(j)
        end do
    end subroutine

    function io_getFileName(default_fname) result(fname)
        character(len=*), optional :: default_fname
        character(len=100) :: fname
        if(.NOT. present(default_fname))default_fname="" 
        write(*,'(3(A))', advance='no') "File (default=", default_fname ,"):"       
        read(*,'(100A)') fname
        if(fname=='') fname = default_fname
        write(*,'(2(A))') "Opening: ", fname
    end function
    
    function io_openFileAsk(default_fname,mode) result(i)
        integer ::i
        Character(len=*) :: mode,default_fname
        do i=7,max_n
            if(nunit(i)) then  
                nunit(i) = .false.
                open(unit=i, file=io_getFileName(default_fname), status=mode) 
                return 
            endif 
        enddo
        i= -1
    end function

    function io_openFile(default_fname,mode) result(i)
        integer ::i
        Character(len=*) :: mode,default_fname
        do i=7,max_n
            if(nunit(i)) then  
                nunit(i) = .false.
                open(unit=i, file=default_fname, status=mode) 
                return 
            endif 
        enddo
        i= -1
    end function
    subroutine io_closeAllFiles
        integer :: i 
        do i=7,max_n
            if(.NOT. nunit(i)) then
                call io_closeFile(i)
            end if
        end do 
    end subroutine
    subroutine io_closeFile(i)
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

    subroutine io_writeXYZ(f,c,xyz)
        integer :: f,i,n,m
        character :: c(:)
        real :: xyz(:)
        n = size(xyz)
        m = size(c)

        call assert(n,m*3,"Dimension error")

        write(f,*) m
        write(f,*) 
        do i = 1,m
            write(f,*) c(i), xyz(i:i+3)
        end do
    end subroutine
    function io_getFileLines(f) result(nlines)
        integer :: f,nlines
        nlines = 0
        do 
            read (f,*,END=10)
            nlines = nlines +1
        end do
        10 rewind(f)
    end function
end module
