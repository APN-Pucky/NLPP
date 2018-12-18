real function eval_pow(self,x)
class(funktion) :: self
real :: x
eval_pow = self%p_f_array(1)%get(x)**self%p_f_array(2)%get(x)
end function
real function eval_minus(self,x)
class(funktion) :: self
real :: x
eval_minus = self%p_f_array(1)%get(x)-self%p_f_array(2)%get(x)
end function
real function eval_plus(self,x)
class(funktion) :: self
real :: x
eval_plus = self%p_f_array(1)%get(x)+self%p_f_array(2)%get(x)
end function
real function eval_mult(self,x)
class(funktion) :: self
real :: x
eval_mult = self%p_f_array(1)%get(x)*self%p_f_array(2)%get(x)
end function
real function eval_div(self,x)
class(funktion) :: self
real :: x
eval_div = self%p_f_array(1)%get(x)/self%p_f_array(2)%get(x)
end function
function base_pow(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(c,2)
C%p_f_array = (/B,A/)
C%p_eval => eval_pow
end function
function base_minus(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(c,2)
C%p_f_array = (/B,A/)
C%p_eval => eval_minus
end function
function base_plus(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(c,2)
C%p_f_array = (/B,A/)
C%p_eval => eval_plus
end function
function base_mult(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
CALL alloc(C,2)
C%p_f_array = (/B,A/)
C%p_eval => eval_mult
end function
function base_div(A,B) result(C)
class(funktion) :: A
type(funktion) :: B,C
call alloc(c,2)
C%p_f_array = (/B,A/)
C%p_eval => eval_div
end function
function plus_f(A,B) result(C)
class(funktion),INTENT(IN) :: A
type(funktion), INTENT(IN) :: B 
type(funktion) :: C
C = B%p_plus(A)
end function
function plus_real_r(A,B) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = A+(d*B)
end function
function plus_real_l(B,A) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = (d*B)+A
end function
function minus_f(A,B) result(C)
class(funktion),INTENT(IN) :: A
type(funktion), INTENT(IN) :: B 
type(funktion) :: C
C = B%p_minus(A)
end function
function minus_real_r(A,B) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = A-(d*B)
end function
function minus_real_l(B,A) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = (d*B)-A
end function
function pow_f(A,B) result(C)
class(funktion),INTENT(IN) :: A
type(funktion), INTENT(IN) :: B 
type(funktion) :: C
C = B%p_pow(A)
end function
function pow_real_r(A,B) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = A**(d*B)
end function
function pow_real_l(B,A) result(C)
class(funktion),INTENT(IN) :: A
real,INTENT(IN) :: B 
type(funktion) :: C
C = (d*B)**A
end function
function mult_f(A,B) result(C)
class(funktion),INTENT(IN) :: A
type(funktion), INTENT(IN) :: B 
type(funktion) :: C
C = B%p_mult(A)
end function
function mult_real_r(A,B) result(C)
class(funktion),INTENT(IN) :: A
real, INTENT(IN) :: B
type(funktion) :: C
C = A
C%scal = C%scal * B
end function
function mult_real_l(B,A) result(C)
class(funktion),INTENT(IN) :: A
real, INTENT(IN) :: B
type(funktion) :: C
C = A
C%scal = C%scal * B
end function
function div_f(A,B) result(C)
class(funktion),INTENT(IN) :: A
type(funktion), INTENT(IN) :: B 
type(funktion) :: C
C = B%p_div(A)
end function

function div_real_r(A,B) result(C)
class(funktion),INTENT(IN) :: A
real, INTENT(IN) :: B
type(funktion) :: C
C = A
C%scal = C%scal/B
end function
function div_real_l(B,A) result(C)
class(funktion),INTENT(IN) :: A
real, INTENT(IN) :: B
type(funktion) :: C
C = d/A*B
end function
