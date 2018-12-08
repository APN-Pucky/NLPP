program main
    use io
    implicit none
    call io_selectLoop((/&
        Select("a) Einfach Euler HO", eeho),&
        Select("b) Korrigierter Euler HO", keho),&
        Select("c) Iterativer Euler HO", ieho),&
        Select("d) Kepler", kepler)&
    /))

    contains 

    subroutine eeho()
        integer :: i,f
        ! vx = z && vz = -x
        real :: dt=0.001,t=100
        real :: x=2.,vx=0.
        !real :: vz=-2,z=0.
        f = io_openFile('data/eeho.dat',"replace")
        do i=1,int(t/dt)
            vx = vx-x *dt
            !vx = z
            x = x+vx *dt
            !vz = -x
            write(f,*) i*dt,x
        end do
        call io_close(f)
    end subroutine

    subroutine keho()
        integer :: i,f
        ! vx = z && vz = -x
        real :: dt=0.001,t=100
        real :: x=2.,vx=0.
        !real :: vz=-2,z=0.
        real :: tvx, tx
        f = io_openFile('data/keho.dat',"replace")
        do i=1,int(t/dt)
            !vx
            tvx = vx-x *dt
            tx = x+tvx *dt
            vx = vx - (x+tx)/2*dt
            print *, x , tx
            !x
            tx = x + vx*dt
            tvx = vx -tx *dt
            x = x + (vx+tvx)/2*dt
            write(f,*) i*dt,x
        end do
        call io_close(f)
    end subroutine

    subroutine ieho()
        integer :: i,j,f,n =2
        ! vx = z && vz = -x
        real :: dt=0.001,t=100
        real :: x=2.,vx=0.
        !real :: vz=-2,z=0.
        real :: tvx, tx
        f = io_openFile('data/ieho.dat',"replace")
        do i=1,int(t/dt)
            !vx
            tx = x
            do j=1,n
                tvx = vx-(x+tx)/2*dt
                tx = x+tvx *dt
            end do
            vx = vx - (x+tx)/2*dt
            print *, x , tx
            !x
            tvx = vx
            do j =1,n
                tx = x + (vx+tvx)/2*dt
                tvx = vx -tx *dt
            end do
            x = x + (vx+tvx)/2*dt
            write(f,*) i*dt,x
        end do
        call io_close(f)
    end subroutine

    subroutine kepler()
        integer :: i,j,f,n,ff
        real :: dt = 0.01,t=100
        real :: x=0.5,vx=0,y=0,vy=1.63,z=0
        real :: tx,tvx,ty,tvy
        f = io_openFile('data/kepler.xyz',"replace")
        ff = io_openFile('data/kepler.dat',"replace")
        do i=1,int(t/dt)
            !vx
            tvx = vx-x/(sqrt((x*x+y*y))**3)*dt
            tx = x+tvx *dt
            vx = vx - (x+tx)/2/(sqrt((x*x+y*y))**3)*dt
            !print *, x , tx

            !vy
            tvy = vy-y/(sqrt((x*x+y*y))**3)*dt
            ty = y+tvy *dt
            vy = vy - (y+ty)/2/(sqrt((x*x+y*y))**3)*dt
 
            !x
            tx = x + vx*dt
            tvx = vx -tx/(sqrt((tx*tx+y*y))**3) *dt

            !y
            ty = y + vy*dt
            tvy = vy -ty/(sqrt((x*x+ty*ty))**3) *dt

            !xy
            x = x + (vx+tvx)/2*dt
            y = y + (vy+tvy)/2*dt
            write(ff,*) i*dt,x,y
            write(f,*) '2'
            write(f,*)
            write(f,*) 'H',x,y,z
            write(f,*) 'Au',z,z,z
        end do
        call io_close(f)
        call io_close(ff)
   end subroutine 
end program
