program main
    use constants
    use eigen
    implicit none
    integer,parameter ::n= 6 ! dimension
    integer :: j=1,i=1
    real, dimension(n,n) :: A! = reshape ( (/1,1,0,1,2,1,0,1,1/), (/3,3/),order=(/2,1/))
    real, dimension(n,n) :: SA! = reshape ( (/1,1,0,1,2,1,0,1,1/), (/3,3/),order=(/2,1/))
    real, dimension(n,n) :: X!= reshape( (/1,0,0,0,1,0,0,0,1/), (/3,3/),order=(/2,1/))
    real, dimension(n,n) :: T=0;
    call hueckel(A)
    call hueckel(SA)
    ! zyklisch?
    A(n,1)=1
    SA(n,1)=1
    A(1,n)=1
    SA(1,n)=1
    ! zyklisch?
    call eigencalc(A,X)
    print * , "A"
    call pr(SA)
    print *, "NA"
    call pr(A)
    print *, "NX"
    call pr(X)
    print *, "----------"
    print *, "check"
    print *, "----------"
    print * , "A*NX ?= NA*NX"
    print * , "A*NX"
    T = matmul(SA,X)
    call pr(T)
    print * , "NX*NA"
    T = matmul(X,A)
    call pr(T)
    print *, "----------"
    print *, "result"
    print *, "----------"

    print *, "Eigenvektoren"
    print *, ((X(i,j), "|",j=1,n),nl,i=1,n)
    print *, "Eigenwerte"
    print *, (A(i,i) ,"|" ,i=1,n)

   
    contains 
    subroutine hueckel(A)
        real :: A(:,:)
        integer :: n,i,j
        n = size(a,1)
        do i =1,n
            do j =1,n
                
               if(i==j+1 .OR. i==j-1) then
                    A(i,j)=1
                else
                    A(i,j)=0
                end if
            end do
        end do
    end subroutine
                
!    do i=1,n
!        L(i,i) = 1
!    end do
!
!    do j=1,n
!        do i=1,n
!            if(i<=j) then
!                t=0
!                do k=1,i-1
!                    t = t+L(i,k)*U(k,j) 
!                end do
!                U(i,j) = A(i,j)-t
!            else
!                t=0
!                do k=1,j-1
!                    t = t+L(i,k)*U(k,j) 
!                end do
!                L(i,j) = (A(i,j)-t)/U(j,j)
!            end if
!        end do
!    end do
!
!    print *, "U"
!    print *, (U(k,:),nl,k=1,n)
!    print *, "L"
!    print *, (L(k,:),nl,k=1,n)
!    LU = matmul(L,U)
!    print *, "LU"
!    print *, (LU(k,:),nl,k=1,n)
!
!    do i =1,n
!        t=0
!        do j=1,i-1 
!            t=t+L(i,j)*y(j)
!        end do
!        y(i) = (b(i) - t)/L(i,i)
!    end do 
!
!    do i =n,1,-1
!        t=0
!        do j=i+1,n
!            t=t+U(i,j)*x(j)
!        end do
!        x(i) = (y(i) - t)/U(i,i)
!    end do 
!
!    print *, "x"
!    print *, x
!
!    t=1        
!    do i=1,n
!        t=t*U(i,i)
!    end do
!    print *,"det(A)"
!    print *, t
    

end program
