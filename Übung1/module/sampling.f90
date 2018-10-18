module sampling
    use mfunktion
    implicit none

    contains 

    subroutine sampling_main()
        implicit none
        type(funktion) :: f,ex,w,xu
        integer,parameter :: n_max_pot = 8
        real :: b=1,a=0.1, alpha = 0.8497369
        real :: max_y = 2.16443469
        integer :: n=100000, i =0
        double precision :: sum =0, int=0, eval
        real :: teval
        real :: p,y
        integer :: n_tot=0,n_hit=0,n_it=1,j

        call init !setup x,dx,d!
        f = x**(-1/3.)+x/10.
        do j=1, n_max_pot
            do i=n_tot,n_it-1
                p = rand(a,b)
                y = rand(0.0,max_y)
                n_tot = n_tot +1
                eval = dble(f%get(p))
                if(y < eval) then
                    n_hit = n_hit + 1
                end if
                sum = sum + eval
            end do
            
            !Mittelwert-Monte-Carlo
            int = sum * (b-a)/n_tot
            write(*,*) "log_10(n)=", j, ": Mittelwert        = ", int
    

            !Neumann-Rejection-Monte-Carlo
            int = real(n_hit)/real(n_tot)*max_y*(b-a)
            write(*,*) "log_10(n)=", j ,": Neumann-Rejection = ", int

            n_it = n_it * 10
        end do 
       
        !Importance Sampling 
        w = x**(-1/3.)*alpha
        xu = ((x-1.+3./2.*alpha)/(3./2.*alpha))**(3./2.)
        n_tot =0
        n_it = 1
        n_hit = 0
        sum = 0
        do j=1, n_max_pot
            do i=n_tot,n_it-1
                p = rand(0.0,1.0)
                n_tot = n_tot +1
                teval = xu%get(p)
                eval = dble(f%get(teval))/dble(w%get(teval))
                sum = sum + eval
            end do
            !Mittelwert-Monte-Carlo
            int = sum /n_tot
            write(*,*) "log_10(n)=", j, ": Importance        = ", int
            n_it = n_it * 10
        end do 
       

        !Simpson verfahren
        !N=1000
        !ex = f*dx
        !write(*,*) "Simpson = ", ex%get((/b,a/))
        call f%dealloc()
    end subroutine
    
    function rand(a,b)
        real :: a,b,rand,x
        call random_number(x)
        rand=x*(b-a)+a
    end function
        

end module



