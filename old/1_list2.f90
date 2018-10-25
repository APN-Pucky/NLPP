module mList
    implicit none
    type :: List
        integer :: n=0,cur=1
        real*8,allocatable :: array(:)
        contains 
        generic :: append => append_real, append_real8
        procedure, pass :: alloc
        procedure, pass :: dealloc
        procedure, pass :: append_real
        procedure, pass :: append_real8
        procedure, pass :: toArray
        !procedure, pass :: tail
    end type

    contains

    function toArray(self) result(a)
        class(List) :: self
        real*8, allocatable,dimension(:) :: a 
        integer :: j
        allocate(a(self%cur-1))
        do j=1,self%cur-1
            a(j) = self%array(j)
        end do
    end function

    subroutine alloc(self,i)
        class(List) :: self 
        integer :: i
        self%n = i
        allocate(self%array(i))
    end subroutine

    subroutine dealloc(self)
        class(List) :: self 
        if(allocated(self%array))deallocate(self%array)
    end subroutine



    subroutine append_real8(self,B)
        class(List) :: self
        real*8, intent(in) :: B
        integer :: j
        type(List) :: ret
        if(self%n .LT. self%cur) then
            call ret%alloc(self%n*2)
            ret%cur = self%cur
            do j=1,self%cur-1
                ret%array(j) = self%array(j)
            end do
            !call self%dealloc()
        else
            ret = self
        endif
        ret%array(ret%cur) = B
        ret%cur=ret%cur+1
    end subroutine 


    function append_real(self,B) result(C)
        class(List), intent(in) :: self
        real,intent(in) :: B
        real*8 :: D
        type(List) :: C
        D=B
        C = self%append_real8(D)
    end function


    subroutine list_test()
        real*8 :: t = 2.2
        real*8 :: t2 = 44
        type(List) :: l
        real*8, allocatable :: arr(:)
        l%append(t)
        l = l + t2
        write(*,*) "l", l%array
        l = l + t2
        l = l + t2
        l = l + t2
        l = l + t2
        l = l + t2
        arr = l%toArray()
        write(*,*) arr
        deallocate(arr)
    end subroutine 
end module
