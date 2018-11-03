module random
    use io
    type RNG
        real :: x=0.1, m=566927
        integer :: a=3612,c=5701
        procedure(n),pointer :: next => random_next
    end type
    abstract interface 
            real function n(r)
                        import RNG
                        class(RNG) :: r
            end function
    end interface
    type(RNG) :: r
contains
    subroutine set_fortran_mode()
        r%next => random_fortran_next
    end subroutine
    subroutine set_rng_mode()
        r%next => random_next
    end subroutine
  
    subroutine random_print_next()
        write(*,*) r%next()   
    end subroutine
    subroutine random_save_next()
        integer :: f,n,i
        real :: t, eps = 1e-1
        real,allocatable :: x(:)
        f = io_OpenFile("data/dat.dat","replace")
        n = io_getInteger(10000)
        allocate(x(n))
        do i=1,n
            x(i) = r%next()
            write(f,*) i,x(i)
        end do
        call io_close(f)
        t = SUM(x)/n
        write(*,*) "<x>=", t
        t = SUM(x*x)/n
        write(*,*) "<x²>=", t
        t = sum(x*cshift(x,1))/n
        write(*,*) "C=", t

        do i=1,n
            t = sum((x-cshift(x,i))**2)/n
            if ( abs(t) < eps) then
                write(*,*) "Periode:", i
            end if
        end do
        deallocate(x)
    end subroutine
    real function random_next(self)
        class(RNG) :: self
        self%x = modulo(self%a*self%x+self%c,self%m)
        random_next = self%x/self%m
    end function
    real function random_fortran_next(self)
        class(RNG) :: self
        call random_number(random_fortran_next)
    end function
       
end module
