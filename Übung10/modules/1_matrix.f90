module matrix
    use constants
    use base
    implicit none


    !do i=1,10
    !    print *, i, s(A)
    !    call sweep(A,X)
    !end do
    !print *, "NA"
    !call pr(A)
    !print *, "NX"
    !call pr(X)
    !print *, "----------"
    !print *, "check"
    !print *, "----------"
    !print * , "A*NX ?= NA*NX"
    !print * , "A*NX"
    !T = matmul(SA,X)
    !call pr(T)
    !print * , "NX*NA"
    !T = matmul(X,A)
    !call pr(T)
    !print *, "----------"
    !print *, "result"
    !print *, "----------"

    !print *, "Eigenvektoren"
    !print *, ((X(i,j), "|",j=1,n),nl,i=1,n)
    !print *, "Eigenwerte"
    !print *, (A(i,i) ,"|" ,i=1,n)

   
    contains 

    !TODO check determinat
    subroutine eigencalc(A,X)
        real :: A(:,:),X(size(A,1),size(A,2))
        integer :: n,i
        call assert(size(A,1),size(A,2)) ! quadrativc only
        n = size(A,1)
        CALL set_id(X)
        call pr(X)

        do i=1,10 !10 sweeps
            print *, i, s(A)
            call sweep(A,X)
        end do
    end subroutine 

    subroutine set_id(A)
        real :: A(:,:)
        integer :: n,i,j
        n = size(A,1)
        do i=1,n
            do j=1,n
                if(i==j) then
                    A(i,j) =1
                else
                    A(i,j)= 0
                end if
            end do
        end do
    end subroutine


    real function s(A) !SUM (?trace?)
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
    subroutine wr(f,A)
        real, dimension(:,:) :: A
        integer :: k,f
        write(f,*) (A(k,:),nl,k=1,size(A,2))
    end subroutine

    function lu_determinant(a) result(det)
        integer :: n
        real :: a(:,:)
        real :: det
        real :: u(size(a,1),size(a,2)),l(size(a,1),size(a,2))
        real :: t
        integer :: k=1,i=1,j=1
        n = size(a,1)

        call assert(n,size(a,2)) !quadratic

        u=0
        l=0

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
  
	    det=1        
	    do i=1,n
	        det=det*U(i,i)
	    end do  
    end function

    function lu_factorization(a,b) result(x)
        integer :: n
        real :: a(:,:)
        real :: b(:)
        real :: u(size(a,1),size(a,2)),l(size(a,1),size(a,2))
        real :: x(size(b)), y(size(b))
        real :: t
        integer :: k=1,i=1,j=1

        n = size(a,1)

        call assert(n,size(a,2)) !quadratic
        call assert(n,size(b))

        u=0
        l=0
        x=0
        y=0

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
     

    end function

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
    

end module
