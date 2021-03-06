module verlet
    use io
    !TODO numeric energy diff
    implicit none

    real,parameter :: morse_d=0.493172 
    real,parameter :: morse_a=1.0883376020965 
    real,parameter :: morse_r0=2.2696744602423 


 
    abstract interface 
        function func_xyz_xyz(xyz,n)
            real :: xyz(3,n)
            real :: func_xyz_xyz(3,n)
        end function
        function func_xyz_a(xyz,n)
            real :: xyz(3,n)
            real :: func_xyz_a
        end function
        subroutine  func_step(x,v,f,m,dt,force,n)
            integer :: n
            real, dimension(3,n) :: x,v,f
            real, dimension(n) :: m
            real :: dt
            procedure(func_xyz_xyz),pointer :: force
        end subroutine
    end interface

    type :: Step
        Character (len=20) :: name
        procedure(func_step),pointer,nopass :: exec
    end type


    contains

    real function total_energy(x,v,m,pot,n)
        integer :: n
        real, dimension(3,n) :: x,v
        real :: m(n)
        procedure(func_xyz_a),pointer :: pot
        total_energy = pot(x,n)+ kinetic_energy(v,m,n)
    end function

    real function kinetic_energy(v,m,n) result(kin)
        integer :: n,i
        real :: m(n)
        real, dimension(3,n) :: v
        Do i=1,n
            kin = kin + m(i)*SUM(v(:,i)**2)
        end do
        kin = kin/2
    end function
        

    subroutine simple_euler_step(x,v,f,m,dt,force,n)
        integer :: n,i
        real, dimension(3,n) :: x,v,f
        real, dimension(n) :: m
        real :: dt
        procedure(func_xyz_xyz),pointer :: force
        x=x+v*dt
        do i=1,n
            v(:,i)=v(:,i)+f(:,i)*dt/(m(i))
        end do

        f = force(x,n)
    end subroutine

    subroutine corrected_euler_step(x,v,f,m,dt,force,n)
        integer :: n,i
        real, dimension(3,n) :: x,v,f,tx,tv,tf,ttx
        real, dimension(n) :: m
        real :: dt
        procedure(func_xyz_xyz),pointer :: force
        !x
        tx=x+v*dt
        tf = force(tx,n)
        do i=1,n
            tv(:,i)=v(:,i)+tf(:,i)*dt/(m(i))
        end do
        tx = x
        x = x+ (v+tv)/2*dt
        !v
        !no change to x (tx is old x)
        do i=1,n
            tv(:,i)=v(:,i)+f(:,i)*dt/(m(i))
        end do
        ttx = tx+tv*dt
        tf = force(ttx,n)
        do i=1,n
            v(:,i)=v(:,i)+(f(:,i)+tf(:,i))/2*dt/(m(i))
        end do

        f = force(x,n)
    end subroutine




    subroutine verlet_step(x,v,f,m,dt,force,n)
        integer :: n,i
        real, dimension(3,n) :: x,v,f
        real, dimension(n) :: m
        real :: dt
        procedure(func_xyz_xyz),pointer :: force

        do i=1,n
            v(:,i)=v(:,i)+f(:,i)*dt/(2*m(i))
        end do
        x=x+v*dt
        f=force(x,n)
        do i=1,n
            v(:,i)=v(:,i)+f(:,i)*dt/(2*m(i))
        end do
    end subroutine

    ! CO Test/Sim

    subroutine test_CO
        integer :: dt(3) = (/10,40,100/)
        real :: r0(3) = (/ 2.5,2.3,3.0/)
        type(Step) :: steps(3) 
        integer i,j,k
        steps = (/ Step("verlet",verlet_step), Step("corrected_euler",corrected_euler_step), &
        Step("simple_euler",simple_euler_step) /)
        do i=1,size(dt)
            do j=1,size(r0)
                do k=1,size(steps)
                    call sim_CO_cms(steps(k),dt(i),r0(j))
                    call sim_CO_lab(steps(k),dt(i),r0(j))
                end do 
            end do 
        end do
    end subroutine

    subroutine sim_CO_cms(s,dt,r0)
        integer,parameter :: n=1
        real :: r0
        integer :: dt
        real :: o_xyz(3,n)!reshape((/ r0,0.,0./),(/ 3,n /))
        real :: c_xyz(3,n)!reshape((/ 0.,0.,0./),(/ 3,n /))
        real :: v(3,n) !reshape((/ 0.,0.,0./),(/ 3,n /))
        real :: f(3,n) !reshape((/0.,0.,0./),(/ 3,n /))
        real :: m(n) = (/12506/) !m_C*m_O/(m_C+m_O)/m_e *[m_e]
        integer :: i,ff,fd
        integer :: max_dt = 10000
        procedure(func_xyz_xyz),pointer :: force => grad_morse_potential_cms
        procedure(func_xyz_a),pointer :: pot => morse_potential_cms
        type(Step) :: s 
        o_xyz = 0.
        c_xyz = 0.
        v = 0.
        f = 0.
        o_xyz(1,1) = r0
        ff = io_openFile("data/CO_cms_" // trim(s%name)//"_dt" //trim(itoa(dt))// "_r" //trim(rtoa(r0))//".xyz","replace")
        fd = io_openFile("data/CO_cms_" // trim(s%name)//"_dt" //trim(itoa(dt))// "_r" //trim(rtoa(r0))//".dat","replace")
        f=force(o_xyz,n)
        do i=1,max_dt,dt
            CALL s%exec(o_xyz,v,f,m,real(dt),force,n)
            CALL io_writeXYZ(ff,(/'O','C'/), (/ o_xyz(:,:), c_xyz(:,:) /))
            write(fd,*) i, sqrt(sum(o_xyz(:,1)**2)), total_energy(o_xyz,v,m,pot,n), v(1,1)
        end do
        call io_closeFile(ff)
        call io_closeFile(fd)
    end subroutine

    subroutine sim_CO_lab(s,dt,r0)
        integer,parameter :: n=2
        real :: r0
        integer :: dt
        real :: oc_xyz(3,n)!reshape((/ 2.5,0.,0.,0.,0.,0./),(/ 3,n /))
        real :: v(3,n) !reshape((/ 0.,0.,0.,0.,0.,0./),(/ 3,n /))
        real :: f(3,n) !reshape((/0.,0.,0.,0.,0.,0./),(/ 3,n /))
        real :: m(n) = (/29164,21895/) !m_O,m_C
        integer :: i,ff,fd
        integer :: max_dt = 10000
        procedure(func_xyz_xyz),pointer :: force => grad_morse_potential_lab
        procedure(func_xyz_a),pointer :: pot => morse_potential_lab
        type(Step) :: s
        oc_xyz = 0.
        v = 0.
        f = 0.
        oc_xyz(1,1)= r0
        ff = io_openFile("data/CO_lab_" // trim(s%name)//"_dt" //trim(itoa(dt))// "_r" //trim(rtoa(r0))//".xyz","replace")
        fd = io_openFile("data/CO_lab_" // trim(s%name)//"_dt" //trim(itoa(dt))// "_r" //trim(rtoa(r0))//".dat","replace")
        f=force(oc_xyz,n)
        do i=1,max_dt,dt
            CALL s%exec(oc_xyz,v,f,m,real(dt),force,n)
            CALL io_writeXYZ(ff,(/'O','C'/), (/ oc_xyz(:,:) /))
            write(fd,*) i, sqrt(sum((oc_xyz(:,1)-oc_xyz(:,2))**2)), total_energy(oc_xyz,v,m,pot,n)
        end do
        call io_closeFile(ff)
        call io_closeFile(fd)
    end subroutine

    ! Morse Potential

    real function morse_potential_cms(xyz,n)
        integer ::n 
        real :: xyz(3,n)
        real :: r,e
        call assert(n,1)

        r= sqrt(sum(xyz(:,1)**2))
        e = exp(-morse_a*(r-morse_r0))
        morse_potential_cms = morse_d*(1-e)**2
    end function
        
    function grad_morse_potential_cms(xyz,n)
        integer :: n
        real :: xyz(3,n)
        real :: grad_morse_potential_cms(3,n)
        real :: e,r

        CALL assert(n,1)

        r= sqrt(sum(xyz(:,1)**2))
        e = exp(-morse_a*(r-morse_r0))
        grad_morse_potential_cms(:,1) = -morse_d*2*morse_a*xyz(:,1)*e*(1-e)/r
    end function


    real function morse_potential_lab(xyz,n)
        integer ::n 
        real :: xyz(3,n)
        real :: r,e
        call assert(n,2)

        r= sqrt(sum((xyz(:,1)-xyz(:,2))**2))
        e = exp(-morse_a*(r-morse_r0))
        morse_potential_lab = morse_d*(1-e)**2
    end function

    function grad_morse_potential_lab(xyz,n)
        integer :: n
        real :: xyz(3,n)
        real :: grad_morse_potential_lab(3,n)
        real :: e,r

        CALL assert(n,2)

        r= sqrt(sum((xyz(:,1)-xyz(:,2))**2))
        e = exp(-morse_a*(r-morse_r0))
        grad_morse_potential_lab(:,1) = -morse_d*2*morse_a*(xyz(:,1)-xyz(:,2))*e*(1-e)/r
        grad_morse_potential_lab(:,2) = morse_d*2*morse_a*(xyz(:,1)-xyz(:,2))*e*(1-e)/r
    end function

end module
