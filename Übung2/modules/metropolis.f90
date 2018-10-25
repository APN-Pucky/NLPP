module metropolis
    use mIO
    use mElement
    use mConstants
    use mMath
    implicit none
    integer :: n
    type(Element),dimension(:),allocatable :: start_l
    type(Element),dimension(:),allocatable :: l

    !parameter
    integer ::nstep=300000
    real*8 :: t=2000,delta=0.12
    contains 
    subroutine metro()
        call io_selectLoop((/&
            Select("Load/Init Data",load_data),&
            Select("Print Data", print_data),&
            Select("Calc & Save Data",calc_data)&
        /))

    end subroutine

    subroutine mcmove(x,y,z,n,moved,e)
        real*8, dimension(n) :: x,y,z,xnew,ynew,znew
        integer :: n,i 
        logical :: moved
        real*8 :: eold,rand,boltz
        real*8, optional :: e

        ! calculate energy at old positions
        if(.NOT. present(e))CALL energy(x,y,z,n,e)       
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
           DO i=1,n
        
              x(i)=xnew(i)
              y(i)=ynew(i)
              z(i)=znew(i)
        
           ENDDO
                 
        ! count accepted trial moves
           !nacc=nacc+1
           moved = .true.
        ELSE
                 
           e=eold
           moved=.false.
                 
        ENDIF
              
        RETURN
        
    end subroutine

    subroutine calc_data()
        integer :: nacc,istep=1,f,i
        logical :: moved
        real*8 :: e=0,c
        real*8,allocatable :: energies(:) 
        nacc=0
        l = start_l
        f= io_openFile("data/au13-out.xyz",'replace') 
        write(f,*) n
        ! nstep MC steps 
        DO istep=1,nstep 

            ! randomly displace atoms
            CALL mcmove(l%x,l%y,l%z,n,moved,e) 
            if(moved) then
                nacc = nacc+1
                write(f,*) e,nacc
                write(f,*) (l(i),nl,i=1,n),n
            endif
        ENDDO
        rewind(f)
        allocate(energies(nacc))
        Do istep=1, nacc
            call io_skip(f,1)
            read(f,*) energies(istep), i
            call io_skip(f,n)
        end do
        c=variance(energies)/kb/t**2
        deallocate(energies)
        call io_close(f)
        write(*,'(A,I0,A,F5.2,A)') "Accepted steps: " , nacc," (", real(nacc)/real(nstep)*100,"%)"
        write(*,*) "c=", c
    end subroutine


    subroutine load_data()
        integer :: b,i,f
        !call load_data("",ttt)
        if(allocated(start_l))deallocate(start_l)
        if(allocated(l))deallocate(l)
        f = io_openFile("data/au13-ini.xyz",'old')
        read(f,*) n,b
        allocate(start_l(n))
        read(f,*) (start_l(i),i=1,n)
        call io_close(f)
    end subroutine

    subroutine print_data()
        integer :: i
        !call load_data("",ttt)
        write(*,*) (l(i),nl,i=1,n)
    end subroutine

end module
