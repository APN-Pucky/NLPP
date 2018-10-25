module list
    implicit none
    type :: Node
        real*8 :: element
        logical :: s = .false.
        type(Node), pointer :: next => null()
        contains 
        generic :: append => append_real, append_real8
        procedure, pass :: hasNext
        procedure, pass :: set
        procedure, pass :: dealloc
        procedure, pass :: remaining
        procedure, pass :: append_real
        procedure, pass :: append_real8
        procedure, pass :: toArray
        !procedure, pass :: tail
    end type
    interface operator(+)
        module procedure append_real8, append_real
    end interface

    contains

    function toArray(self) result(a)
        class(Node) :: self
        real*8, allocatable,dimension(:) :: a 
        integer :: n,i
        type(Node) :: cur
        
        write(*,*) "alloc"

        n = 1+self%remaining()
        allocate(a(n))
        cur = self
        Do i=1,n
            a(i) = cur%element
            cur = cur%next
        end do
    end function

    logical function hasNext(self)
        class(Node) :: self
        hasNext = associated(self%next)
    end function

    integer recursive function remaining(self) result(n)
        class(Node) :: self
        n= 0
        if(self%hasNext()) then
            write(*,*) self%element
            n= 1+self%next%remaining()
        end if
    end function  

    recursive function append_real8(self,B) result(d)
        class(Node) :: self
        real*8, intent(in) :: B
        type(Node),allocatable,target :: C
        type(Node) :: d,f
        if(self%hasNext()) then 
            f=self%next%append_real8(B)
        else
            allocate(c)
            call C%set(B)
            self%next => C
            write(*,*) self%element, "=>", C%element
        end if
        d%next = self%next
        d%s  = self%s
        d%element = self%element
        write(*,*) "d" , d%next%element
        !if(self%s)C%next => self
    end function 


    function append_real(self,B) result(C)
        class(Node),target, intent(in) :: self
        real,intent(in) :: B
        real*8 :: D
        type(Node) :: C
        D=B
        C = self%append_real8(D)
    end function

    subroutine set(self,v)
        class(Node) :: self
        real*8 :: v
        self%element = v
        self%s = .true.
    end subroutine

    subroutine dealloc(self)
        class(Node) :: self
        if(self%hasNext()) then
            call self%next%dealloc()
            deallocate(self%next)
        endif
    end subroutine

    subroutine list_test()
        real*8 :: t = 2.2
        real*8 :: t2 = 44
        type(Node) :: list
        real*8, allocatable :: arr(:)
        write(*,*) list%hasNext()
        list%element = 7.
        list =list%append(t)
        write(*,*) "l", list%next%element
        list = list + t2
        write(*,*) "l", list%element
        list = list + t2
        arr = list%toArray()
        write(*,*) arr
        deallocate(arr)
    end subroutine 
end module
