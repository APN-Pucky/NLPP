program main
    use constants
    implicit none

    real, dimension(3,3) :: A = reshape ( (/1,1,0,1,2,1,0,1,1/), (/3,3/),order=(/2,1/))
    real, dimension(3,3) :: SA = reshape ( (/1,1,0,1,2,1,0,1,1/), (/3,3/),order=(/2,1/))
    real, dimension(3,3) :: X= reshape( (/1,0,0,0,1,0,0,0,1/), (/3,3/),order=(/2,1/))
    real, dimension(3,3) :: T=0;
    integer :: j=1,i=1,n=3

    do i=1,10
        print *, i, s(A)
        call sweep(A,X)
    end do
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

    real function s(A)
        real :: A(:,:)
        integer :: n,i,j
        s = 0
        n = size(A,1)
        do i=1,n
            do j=1,n    
                if(i/=j)s=s + A(i,j)*A(i,j)
            end do
        end do
    end function

    subroutine sweep(A,X)
        real :: A(:,:)
        real :: X(size(A,1),size(A,2))
        integer :: n,i,j
        n = size(A,1)
    
        !call pr(A)
        do i=1,n
           do j =i,n
                if(i/=j)call rot(A,X,i,j)
           end do
        end do
    end subroutine

    subroutine rot(A,X,p,q)
        real :: A(:,:)
        real :: X(size(A,1),size(A,2))
        real :: NX(size(A,1),size(A,2))
        real :: NA(size(A,1),size(A,2))
        integer :: p,q
        integer :: n,i,j
        real :: theta 
        real :: t 
        real :: c
        real :: s 
        n = size(A,1)
        theta = (A(q,q)-A(p,p))/(2*A(p,q))
        t =sign(1.,theta)/(abs(theta)+sqrt(theta*theta+1))
        c = 1/sqrt(1+t*t)
        s= c*t
        NA = A
        NX = X
        do i=1,n
            do j=1,n
                if(j==p) then
                    if(i==p) then 
                        NA(p,p)= a(p,p)-t*A(p,q)
                    else if(i==q) then
                        NA(p,q)=0
                        NA(q,p)=0
                    else 
                        NA(i,j) = c*A(i,p)-s*A(i,q)
                        NA(j,i) = NA(i,j)
                    end if
                    NX(i,p) = c*X(i,p)-s*X(i,q)
                else if(j==q) then
                    if(i==p) then 
                        NA(p,q)=0
                        NA(q,p)=0
                    else if (i==q) then
                        NA(i,j)= A(q,q)+t*A(p,q)
                    else
                        NA(i,j) = c*A(i,q)+s*A(i,p)
                        NA(j,i) = NA(i,j)
                    end if
                    NX(i,q) = c*X(i,q)+s*X(i,p)
                end if
            end do
        end do
    
        A=NA
        X=NX
    end subroutine

    !Array print
    subroutine pr(A)
        real, dimension(:,:) :: A
        integer :: k
        print *, (A(k,:),nl,k=1,size(A,2))
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
