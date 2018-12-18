module metal
    use mfunktion
    use matrix
    implicit none
    contains    

    function ijtol(i,j,n,m) result(l)
        integer :: i,j,l,n,m
        l = i+ (m-1-j)*(n-2)-1
    end function

    function ltoi(l,n,m) result(i)
        integer :: i,l,n,m
        i = modulo(l-1,n-2)+2
    end function 

    function ltoj(l,n,m) result(j)
        integer :: j,l,n,m
        j = -(l+1-ltoi(l,n,m))/(n-2)+m-1
    end function

    function sim_metal(length,width,n,m,borders) result(grid)
        integer :: n,m
        integer :: k
        real :: length,width
        type(funktion) :: borders(4)
        real,dimension(n,m) :: grid
        real,dimension((n-2)*(m-2),(n-2)*(m-2)) :: lgs 
        real,dimension((n-2)*(m-2)) :: vec,res
        integer :: r,c,i,j,ti,tj

        k=(n-2)*(m-2)

        !fill lgs and vec
        lgs = 0
        vec = 0
        do r=1,k
            i = ltoi(r,n,m)
            j = ltoj(r,n,m)
            call assert(i/=0)
            call assert(i/=1)
            call assert(j/=0)
            call assert(j/=1)
            call assert(i/=n)
            call assert(j/=m)
            do c=1,k
                ti = ltoi(c,n,m)
                tj = ltoj(c,n,m)
                if (ti==i .and. tj==j) then
                    lgs(r,c)=4
                    if(i==2) then
                        vec(r) = vec(r) + eval_border(i-1,j,length,width,n,m,borders)
                    end if
                    if(j==2) then
                        vec(r) = vec(r) + eval_border(i,j-1,length,width,n,m,borders)
                    end if
                    if(i==n-1) then 
                        vec(r) = vec(r) +  eval_border(i+1,j,length,width,n,m,borders)
                    end if
                    if(j==m-1) then
                        vec(r) = vec(r) + eval_border(i,j+1,length,width,n,m,borders)
                    end if
                else if (abs(ti-i)==1 .and. abs(tj-j)==0) then 
                    lgs(r,c)=-1
                else if (abs(tj-j)==1 .and. abs(ti-i)==0) then 
                    lgs(r,c)=-1
                end if
            end do
        end do

        res = lu_factorization(lgs,vec)
        print *, res
        print *, res(ijtol(2,2,n,m))
        do i=1,n
            do j=1,m
                    if(i==1 .OR. j==1 .OR. i == n .OR. j==m) then
                        grid(i,j) = eval_border(i,j,length,width,n,m,borders)
                    else
                        grid(i,j) = res(ijtol(i,j,n,m))
                    end if
            end do
        end do
                
        call pr(grid)

    end function

    function eval_border(i,j,length,width,n,m,borders) result(r)
        integer :: n,m
        integer :: i,j
        real :: length,width
        type(funktion) :: borders(4)
        real ::r 

            r = 0.
            if(i==1) then
                 r = r + borders(1)%get(length*(j-1)/(n-1))
            else if(j==1) then
                r = r + borders(2)%get(width*(i-1)/(m-1))
            else if(i==n) then 
                 r = r + borders(3)%get(length*(j-1)/(n-1))
            else if(j==m) then
                 r = r + borders(4)%get(width*(i-1)/(m-1))
            end if
    end function

    subroutine test_metal
        integer,parameter :: n=5,m=5
        real :: grid(n,m)=0
        real :: length=0.5,width=0.5
        type(funktion) :: borders(4)
        borders(1) = nil
        borders(2) = nil
        borders(3) = x*200
        borders(4) = x*200
        grid= sim_metal(length,width,n,m,borders)
    end subroutine

end module
