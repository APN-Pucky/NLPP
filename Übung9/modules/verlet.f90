module verlet
    use io
    !TODO numeric energy diff
    implicit none
    abstract interface 
        function func_xyz(xyz,n)
            real :: xyz(3,n)
            real :: func_xyz(3,n)
        end function
    end interface

    contains

    subroutine verlet_step(x,v,f,m,dt,force)
        real, dimension(3,n) :: x,v,f
        real :: m,dt
        procedure(func_xyz),pointer :: force
        v=v+f*dt/(2*m)
        x=x+v*dt
        f=force(x)
        v=v+f*dt/(2*m)
    end subroutine

    subroutine test_CO
        integer :: n=2
        real :: o_xyz(3)= (/ 2.5,0.,0./)
        real :: c_xyz(3)= (/ 0.,0.,0./)
        real :: v(3) = (/ 0.,0.,0./)
        real :: f(3) = (/0.,0.,0./)
        real :: m = 1. !TODO set this
        integer :: i,ff,dt = 20
        integer :: max_dt = 10000
        procedure(func_xyz),pointer :: force => grad_morse_potential
        ff = io_openFile('data/CO.xyz',"replace")
        f=force(xyz)
        do i=1,max_dt
            CALL verlet_step(xyz,v,f,m,real(dt),force)
            CALL io_writeXYZ((/'C','O'
        end do
        
        
    end subroutine

    function grad_morse_potential(xyz,n)
        real :: xyz(3,n)
        real :: grad_morse_potential(3,n)
        real :: d,a,r0,e,r
        integer :: n,i
        CALL assert(n,2)
        d=0.493172 
        a=1.0883376020965 
        r0=2.2696744602423 

        do i=1,n
            r= sqrt(sum(xyz(:,i)**2))
            e = exp(-a*(r-r0))
            grad_morse_potential(:,i) = -d*2*a*xyz(:,i)*e*(1-e)/r
        end do
    end function
        

end module
