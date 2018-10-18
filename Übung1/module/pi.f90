module pi
    implicit none
    contains 
    subroutine pi_main()
        integer,parameter :: n_max_tot=7 ! 10**n_max_tot iter
        integer :: n_it =1

        integer :: n_hit=0
        integer :: n_tot=0

        real :: npi=0 ! pi value
        real :: x,y
        integer :: i,j
        
        do j=1,n_max_tot 
            do i=n_tot,n_it-1
                call random_number(x)
                call random_number(y)
                n_tot = n_tot +1
                if(x*x+y*y<=1) then
                    n_hit = n_hit +1
                end if
            end do

            npi = real(n_hit)/real(n_tot)*4
            write(*,*) "n_tot = ", n_tot, " pi = ", npi
            n_it = n_it*10
        end do
    end subroutine
end module 



