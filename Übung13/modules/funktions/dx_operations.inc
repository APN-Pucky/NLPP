real function eval_dx_div(self,x) !Nummerische Differnezieren
class(funktion) :: self
real :: x
eval_dx_div = (self%p_f_array(1)%get(x+h) - self%p_f_array(1)%get(x))/h !Anwenden der Vorwaertz-Differenz-Methode
end function

real function eval_dx_mult(self,x) !Nummerisches Integrierer nach Simpson Verfahren
class(funktion) :: self
real,dimension(:) :: x !only 2 param!
real :: h
integer :: i
h = (x(1)-x(2))/(N-1) !Schrittweite
eval_dx_mult = 0
eval_dx_mult = eval_dx_mult + self%p_f_array(1)%get(x(1))
do i=1,N-2,1
eval_dx_mult = eval_dx_mult + (MOD(i,2)+1)*2*self%p_f_array(1)%get(x(2)+h*(i)) 
end do
eval_dx_mult = eval_dx_mult + self%p_f_array(1)%get(x(2))
eval_dx_mult = h*eval_dx_mult/3
end function

function dx_div(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(C,1)
C%p_f_array = (/B/)
C%p_eval => eval_dx_div
end function

function dx_mult(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(C,1)
C%p_f_array = (/B/)
C%p_eval => NULL()
C%p_eval_a => eval_dx_mult
end function

