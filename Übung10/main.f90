program main
    use io
    use metal
    use matrix
    use mfunktion
    implicit none
    call init
    call io_selectLoop((/&
        Select("Test a", a),&
        Select("Test b", b),&
        Select("Test Metal", test_metal)&
    /))

    contains 

    subroutine a()
        integer :: l,n=5,m=5,i=0,j=0
        l = ijtol(i,j,n,m)
        print *,l
        print * , ltoi(l,n,m)
        print * , ltoj(l,n,m)
         print * , ltoi(3,n,m)
        print * , ltoj(3,n,m)

        do i = 1,n
                print *, (ijtol(i,j,n,m),j=1,m)
        end do
       do l=1,(n-2)*(m-2)
            print * , ltoi(l,n,m)
            print * , ltoj(l,n,m)
        end do
     
    end subroutine

    subroutine b()
        type(funktion) :: p
        integer :: i
        p = sinus(2*x)
        do i=1,10
            print *, p%get(real(i))
        end do
    end subroutine
end program
