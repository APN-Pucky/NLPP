program main
    use io
    implicit none
    real,parameter :: pi = 3.14159265359
    real, allocatable :: v(:,:,:), x(:,:,:), auto(:), tauto(:),vert(:),dd(:)
    complex, allocatable :: spauto(:), sptauto(:)
    real :: deltat,diff,dr,rmax,r,rho,rrho,rrrho,dist
    integer :: ierror(4), i ,n,d,steps, f(4),j,o,k,nn,jj
    Character :: c

    !INIT
    n=96
    d=3
    steps=2220
    deltat=20.
    o=n/3
    dr = 0.1
    rmax = 9.8652/2/0.529177
    rrrho = real(o)/((2*rmax)**3)


    f(1) = io_openFile("data/traj.xyz", "old")
    f(2) = io_openFile("data/vert.dat", "replace")

    steps = io_getFileLines(f(1))/(N+2)

    !ALLOC
    allocate(dd(d),vert(nint(rmax/dr)),v(n,d,steps-1),x(n,d,steps),auto(steps-1),tauto(steps-1),spauto(steps-1),sptauto(steps-1))


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

    !calc
    do k=1,nint(rmax/dr)
        r = k*dr
        rrho = 0.
        do i=1, steps
            rho = 0.
            do j=1,o
               nn=0
               do jj=1,o
                        dd = x(j,:,i) - x(jj,:,i)
                        dd = dd - 2*rmax*nint(dd/rmax/2)
                        dist = sqrt(sum(dd**2))
                        if(dist<=r .and. dist>r-dr) then
                            nn = nn+1
                        end if
               end do 
               rho = rho + real(nn)
            end do
            rrho = rrho + rho
        end do
        vert(k) = rrho/steps/(4.*pi/3. *(r**3-(r-dr)**3))/o/rrrho
        write(f(2),*) r, vert(k)
    end do


    
end program
