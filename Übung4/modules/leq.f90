module leq
    use constants
    use sort
contains

    subroutine test_leq()
        integer :: n=3
        real, dimension(3,4) :: A = reshape( (/4.,2.,-1., -9.,-4.,2., 2.,4.,2., 2.,3.,1./),(/3,4/))
        real, dimension(3,4) :: B = reshape( (/3.03,-3.03,6.11, -12.1,12.1,-14.2, 14.,-7.,21., -119.,120.,-139./),(/3,4/))

        print *, "A) no_pivot"
        print *, leq_solve_no_pivot(A)
        print *, "B) pivot"
        print *, leq_solve(B)
    end subroutine


    function leq_solve(m) result(r)
        implicit none
        real,dimension(:,:) :: m ! n x (n+1)
        real, dimension(size(m,1)) :: r
        real :: d
        integer :: i,j,k,n
        n=size(m,1)
        call pivot(m)
        do i=1,n
            do j=1,n
                if (i/=j) then
                    !print *, (m(k,:),nl, k=1,n)
                    m(j,:) = m(j,:) - m(i,:)/m(i,i)*m(j,i) 
                end if
            end do
        end do
        !print *, (m(i,:),nl, i=1,n)
        r = (/(m(i,n+1)/m(i,i),i=1,n)/)
    end function

    subroutine pivot(m)
        real,dimension(:,:) :: m ! n x (n+1)
        real, dimension(size(m,1)) :: t
        integer :: n,i
        integer :: p(size(m,1)) 
        p= (/ (i,i=1,size(m,1)) /) ! permutation
        n=size(m,1)
        t = maxval(ABS(m(:,:n)),2)
        call quicksort(t,p)
        m= reshape((/ (m(p(i),:),i=n,1,-1)/),shape(m),order=(/2,1/))
    end subroutine

    function leq_solve_no_pivot(m) result(r)
        implicit none
        real,dimension(:,:) :: m ! n x (n+1)
        real, dimension(size(m,1)) :: r
        real :: d
        integer :: i,j,k,n
        n=size(m,1)
        do i=1,n
            do j=1,n
                if (i/=j) then
                    !print *, (m(k,:),nl, k=1,n)
                    m(j,:) = m(j,:) - m(i,:)/m(i,i)*m(j,i)
                end if
            end do
        end do
        !print *, (m(i,:),nl, i=1,n)
        r = (/(m(i,n+1)/m(i,i),i=1,n)/)
    end function

end module
