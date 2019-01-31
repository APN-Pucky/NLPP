program main
    use autocorrelation
    use dfourier
    use io
    implicit none
    real, allocatable :: v(:,:,:), x(:,:,:), auto(:), tauto(:)
    complex, allocatable :: spauto(:), sptauto(:)
    real :: deltat,diff,s
    integer :: ierror(4), i ,n,d,steps, f(4),j
    Character :: c
    Character :: cc(6)

    !INIT
    n=96
    d=3
    steps=2220
    deltat=20.

    f(1) = io_openFile("data/traj.xyz", "old")
    f(2) = io_openFile("data/correlation.dat", "replace")
    f(3) = io_openFile("data/spectrum.dat", "replace")
    f(4) = io_openFile("data/verschieb.dat", "replace")

    steps = io_getFileLines(f(1))/(N+2)
    print * , steps

    !ALLOC
    allocate(v(n,d,steps-1),x(n,d,steps),auto(steps-1),tauto(steps-1),spauto(steps-1),sptauto(steps-1))


    !READ
    do i = 1, steps
        read(f(1),*)
        read(f(1),*)
        do j=1,N
            read (f(1), *, iostat=ierror(1)) c, x(j,:,i)
        end do
        if (ierror(1) /=0) then
            print *, "I Error"
            exit
        end if
    end do

    !calc v
    do i=1, steps-1
        do j=1,N
            v(j,:,i) = (x(j,:,i)-x(j,:,i+1))/deltat
        end do
    end do
    !correlate
    call correlate(v,auto)
    call timecorrelate(v,tauto)
    !Output
    do i=1, steps
       write(f(2),*,iostat=ierror(2)) deltat*(i-1), auto(i), tauto(i)
       if (ierror(2)/=0) then
          write(*,*) 'O Error'
          exit
       end if
    end do

    !dft
    call dft((1.0,0.0)*Auto,SpAuto)
    call dft((1.0,0.0)*Tauto,SpTauto)
 
    !Output
    do i=1, steps
       write(f(3),*,iostat=ierror(3)) 1./(steps*(137*0.52918e-8*deltat))*(i-1), abs(SpAuto(i)), abs(SpTauto(i))
       if (ierror(3)/=0) then
          write(*,*) 'O Error'
          exit
       end if
    end do

    !diffusion
    diff = sum(auto)/3*deltat
    print *, "d=", diff
    !diff = sum(tauto)/3*20
    !print *, "d_t=", diff


    !diffusion2
    do i=1, steps
        s=0.
        do j=1,N
            s= s+ sum((x(j,:,i)-x(j,:,1))**2)
        end do
        write(f(4),*,iostat=ierror(4)) deltat*(i-1),s/N
        if (ierror(4)/=0) then
          write(*,*) 'O Error'
          exit
       end if
    end do
        
    
end program
