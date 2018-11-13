module sort
contains
recursive subroutine quicksort(a,b)
  implicit none
  real :: a(:)
  integer,optional :: b(*)
  real x, t
  integer :: first = 1, last
  integer i, j,k

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     if (present(b)) then 
          k = b(i);  b(i) = b(j);  b(j) = k
     end if


     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1),b(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last),b(j + 1 : last))
end subroutine quicksort

end module
