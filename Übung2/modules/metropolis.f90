module metropolis
    use io
    use mElement
    implicit none
    integer n
    type(Element),dimension(:),allocatable :: l
    contains 

    subroutine metro()
        call io_selectLoop((/&
            Select("Load Data",load_data),&
            Select("Print Data", print_data)&
        /))

    end subroutine

    subroutine load_data()
        integer :: b,ierror,i
        !call load_data("",ttt)
        if(allocated(l))deallocate(l)
        open(unit=9, file=io_getFileName("data/au13-ini.xyz"), status='old', IOSTAT=ierror)
        read(9,*) n,b
        allocate(l(n))
        read(9,*) (l(i),i=1,n)
    end subroutine

    subroutine print_data()
        integer :: i
        !call load_data("",ttt)
        write(*,*) (l(i),char(10),i=1,n)
    end subroutine

end module
