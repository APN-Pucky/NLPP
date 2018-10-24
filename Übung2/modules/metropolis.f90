module metropolis
    use io
    use mElement
    use constants
    implicit none
    integer n
    type(Element),dimension(:),allocatable :: l
    contains 
    subroutine metro()
        call io_selectLoop((/&
            Select("Load Data",load_data),&
            Select("Print Data", print_data),&
            Select("Calc Data",calc_data)&
        /))

    end subroutine

    subroutine mcmove(x,y,z,n,nacc,e)
        real*8, dimension(n) :: x,y,z,xnew,ynew,znew
        integer :: n,nacc,istep,i 
        real*8 :: e,eold,rand,t=300,delta=0.5,boltz

        ! calculate energy at old positions
        CALL energy(x,y,z,n,e)       
        eold=e
        
        ! random displacement of atoms: 
        DO i=1,n   
        
           CALL random_number(rand)
           xnew(i)=x(i)+delta*(rand-0.5d0)
           CALL random_number(rand) 
           ynew(i)=y(i)+delta*(rand-0.5d0)
           CALL random_number(rand)	
           znew(i)=z(i)+delta*(rand-0.5d0)
        
        ENDDO
        
        ! calculate energy at new positions
        CALL energy(xnew,ynew,znew,n,e)
        
        ! calculate Boltzmann factor: 
        ! kb = Boltzmann constant, t = temperature
        boltz=EXP(-(e-eold)/(kb*t)) 
                                        
        ! accept or reject new positions
        CALL random_number(rand) 
        IF (rand.LT.boltz) THEN   
           Write(*,*) "hit" 
           DO i=1,n
        
              x(i)=xnew(i)
              y(i)=ynew(i)
              z(i)=znew(i)
        
           ENDDO
                 
        ! count accepted trial moves
           nacc=nacc+1
              
        ELSE
                 
           e=eold
                 
        ENDIF
              
        RETURN
        
    end subroutine

    subroutine calc_data()
        integer :: nacc=0,nstep=1,istep=1
        real*8 :: e=0
        ! nstep MC steps 
        DO istep=1,nstep 

            ! randomly displace atoms
            CALL mcmove(l%x,l%y,l%z,n,nacc,e) 
        ENDDO
    end subroutine


    subroutine load_data()
        integer :: b,ierror,i
        !call load_data("",ttt)
        if(allocated(l))deallocate(l)
        open(unit=9, file=io_getFileName("data/au13-ini.xyz"), status='old', IOSTAT=ierror)
        read(9,*) n,b
        allocate(l(n))
        read(9,*) (l(i),i=1,n)
    end subroutine

    subroutine print_data()
        integer :: i
        !call load_data("",ttt)
        write(*,*) (l(i),char(10),i=1,n)
    end subroutine

end module
