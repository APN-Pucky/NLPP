program main
    use io
    use funktions
    use constants
    use analysis
    implicit none
    call init
    call io_selectLoop((/&
        Select("Test a", test_a),&
        Select("Test DFT", test_dft),&
        Select("Test CO", test_co)&
    /))

    contains 

    subroutine test_dft()
        type(funktion) :: p
        real :: t_start=-1.2
        real :: t_end=1.2
        integer :: i,f
        complex,allocatable :: r(:)

        dt = 0.02
        N =  int((t_end - t_start)/dt)
        allocate(r(N+1))
        p = 1-heaviside(x**2-1)
        r = dft(p,t_start)

        f=io_openFile("data/real.dat","replace")
        write( f,*) ( (i-1)*2*pi/N/dt, real(r(i)*dt/sqrt(2*pi)),nl,i=1,N+1)
        call io_closeFile(f)
        f=io_openFile("data/img.dat","replace")
        write( f,*) ( (i-1)*2*pi/N/dt, aimag(r(i)*dt/sqrt(2*pi)),nl,i=1,N+1)
        call io_closeFile(f)
        f=io_openFile("data/2xfourier.dat","replace")
        r = idft(r,t_start)
        write( f,*) ( i*dt+t_start, real(r(i)),nl,i=1,N)
        call io_closeFile(f)

        deallocate(r)
        call io_closeFile(f)
    end subroutine

    subroutine test_co()
        integer :: f
        integer :: i
        real :: v(4,1000)
        complex,allocatable :: r(:)
        f = io_openFile("data/CO_cms_verlet_dt10_r2.50000000.dat","old")
        read(f,*) v(:,:)
        call io_closeFile(f)
        dt=10
        N=1000
        allocate(r(N+1))
        r=dft(cmplx(v(4,:),0.),1.)

        f=io_openFile("data/co-f-real.dat","replace")
        write( f,*) (i*1./N/dt/(137*5.29e-9), abs(real(r(i)*dt/sqrt(2*pi))),nl,i=1,N+1)
        call io_closeFile(f)

        print *, v(4,:)
        print *, r
        deallocate(r)
    end subroutine


    subroutine test_a()
        type(funktion) :: p
        integer :: i
        p = 1-heaviside(x**2-1)
        do i=1,10
            print *,real(i)/5, p%get(real(i)/5)
        end do
    end subroutine
end program
