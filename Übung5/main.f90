program main
    use constants
    implicit none

    real, dimension(3,3) :: A = reshape ( (/4,-9,2,2,-4,4,-1,2,2/), (/3,3/),order=(/2,1/))
    real, dimension(3,3) :: U = 0
    real, dimension(3,3) :: L = 0
    real, dimension(3,3) :: LU = 0
    real, dimension (3)  :: x=0,y=0, b = (/2,3,1/)
    real :: t
    integer :: k=1,i=1,j=1,n=3

    print * , "A"
    print *, (A(k,:),nl,k=1,n)
    do i=1,n
        L(i,i) = 1
    end do

    do j=1,n
        do i=1,n
            if(i<=j) then
                t=0
                do k=1,i-1
                    t = t+L(i,k)*U(k,j) 
                end do
                U(i,j) = A(i,j)-t
            else
                t=0
                do k=1,j-1
                    t = t+L(i,k)*U(k,j) 
                end do
                L(i,j) = (A(i,j)-t)/U(j,j)
            end if
        end do
    end do

    print *, "U"
    print *, (U(k,:),nl,k=1,n)
    print *, "L"
    print *, (L(k,:),nl,k=1,n)
    LU = matmul(L,U)
    print *, "LU"
    print *, (LU(k,:),nl,k=1,n)

    do i =1,n
        t=0
        do j=1,i-1 
            t=t+L(i,j)*y(j)
        end do
        y(i) = (b(i) - t)/L(i,i)
    end do 

    do i =n,1,-1
        t=0
        do j=i+1,n
            t=t+U(i,j)*x(j)
        end do
        x(i) = (y(i) - t)/U(i,i)
    end do 

    print *, "x"
    print *, x
    
    print *, "y"
    print *, y

    t=1        
    do i=1,n
        t=t*U(i,i)
    end do
    print *,"det(A)"
    print *, t
    

end program
