MODULE mfunktion
        IMPLICIT NONE
        type funktion
                procedure(func), pointer,nopass :: p_f => id
                real :: scal=1
                procedure(rfunc), pointer :: p_eval => base_func
                procedure(rafunc), pointer :: p_eval_a => NULL()
                procedure(fffunc), pointer :: p_div => base_div
                procedure(fffunc), pointer :: p_minus => base_minus
                procedure(fffunc), pointer :: p_plus => base_plus
                procedure(fffunc), pointer :: p_pow => base_pow
                procedure(fffunc), pointer :: p_mult => base_mult
                type(funktion), pointer, dimension(:) :: p_f_array
                logical :: is_allocated = .false.
                CONTAINS
                generic :: get => evaluate,evaluate_a
                procedure, pass :: dealloc
                procedure, pass :: base_func
                procedure, pass :: evaluate
                procedure, pass :: evaluate_a
        end type
        interface nst_newton
                module procedure nst_newton_f,nst_newton_fdf,nst_bisection
        end interface
        abstract interface
                real function func(x)
                        real :: x
                end function 
                real function rfunc(rf, x)
                        import funktion
                        class(funktion) :: rf
                        real :: x
                end function
                real function rafunc(rf, x)
                        import funktion
                        class(funktion) :: rf
                        real,dimension(:) :: x
                end function

                function fffunc(f1,f2)
                        import funktion
                        class(funktion) :: f1
                        type(funktion) :: f2,fffunc
                end function
        end interface

        interface operator(*)
                module procedure mult_f, mult_real_l, mult_real_r
        end interface
        interface operator(/)
                module procedure div_f, div_real_l, div_real_r
        end interface
        interface operator(**)
                module procedure pow_f, pow_real_l, pow_real_r
        end interface
        interface operator(+)
                module procedure plus_f, plus_real_l, plus_real_r
        end interface
        interface operator(-)
                module procedure minus_f, minus_real_l, minus_real_r
        end interface


        type(funktion) :: x,dx,d
        real :: e = exp(1.)
        real :: h = sqrt(1e-7)
        real :: eps = 1e-5
        integer :: N = 5 !odd!
        integer :: max_iter = 10000

        CONTAINS 
                include 'base_operations.inc'
                include 'dx_operations.inc'

                subroutine alloc(self,n)
                        type(funktion) :: self
                        integer :: n
                        allocate(self%p_f_array(n))
                        self%is_allocated = .true.
                end subroutine
                recursive subroutine dealloc(self)
                        class(funktion) :: self
                        integer n,i
                        if(.not. self%is_allocated) return
                        n = size(self%p_f_array)
                        DO i=1,n
                                CALL self%p_f_array(i)%dealloc
                        end do
                        deallocate(self%p_f_array)
                        self%is_allocated = .false.
                end subroutine
                        
                function create_funktion(f)
                        type(funktion) :: create_funktion
                        procedure(func), pointer :: f
                        !create_funktion%p_eval => base_func
                        create_funktion%p_f => f
                end function
                real function base_func(self,x)
                        real :: x
                        class(funktion) :: self
                        base_func = self%p_f(x)
                end function

                real function nst_newton_fdf(f,df,x) ! Newton NST mit gegebener Ableitung!
                        real :: x
                        type(funktion) :: f,df,t
                        integer :: i
                        nst_newton_fdf = x
                        t = f/(df) !Verschiebung des x um delta_x=t(x)! 
                        DO i=1,max_iter
                                if(abs(f%get(nst_newton_fdf))<eps)return !Abbruch, wenn nahe genug an Null!
                                nst_newton_fdf = nst_newton_fdf-t%get(nst_newton_fdf)
                        end do
                        nst_newton_fdf = 0
                        nst_newton_fdf = log(nst_newton_fdf) !ERROR No NST!
                end function

                real function nst_newton_f(f,x) !Newton NST mit numerischer Ableitung!
                        real :: x
                        type(funktion) :: f
                        nst_newton_f = nst_newton(f,f/dx,x) !nummerische Ableitung wird als anayltische verwendet!
                end function

                real function nst_bisection(f,x1,x2)
                        real :: x1,x2,bx1,bx2,fx1,fx2,t=0.,ft
                        real :: fm
                        type(funktion) :: f
                        integer :: i
                        bx1 = x1 !avoid pointer/memory issues!
                        bx2 = x2
                        nst_bisection = log(t)
                        fx1 = f%get(bx1)
                        fx2 = f%get(bx2)
                        DO i =1,max_iter
                                if(sign(1.,fx1) .EQ. sign(1.,fx2)) return !keine Gewissheit für NST!
                                t = (bx1+bx2)/2.
                                ft = f%get(t)
                                if(abs(bx1-t) < eps) return !main exit delta x < eps!
                                if(sign(1.,ft) .EQ. sign(1.,fx1)) then
                                        bx1=t
                                        fx1=ft
                                        nst_bisection = t
                                        if(abs(fx1) < eps) return !delta f(x) < eps!
                                else
                                        bx2=t
                                        fx2=ft
                                        nst_bisection = t
                                        if(abs(fx2) < eps) return !delta f(x) < eps!
                                end if

                        end do
                        t = 0
                        nst_bisection = log(t) !ERROR OR max_iter too low!
                end function

                real function minimum(f,x1,x2)
                        !Bisection => local min
                        real :: x1,x2,bx1,bx2, fx1,fx2,t=0.
                        real :: fm
                        type(funktion) :: f
                        integer :: i
                        bx1 = x1 !avoid pointer issues
                        bx2 = x2
                        fm = -log(t)
                        DO i=1,max_iter
                                fx1 = f%get(bx1)
                                fx2 = f%get(bx2)
                                if(fx1<fm) then 
                                        minimum = bx1
                                        fm = fx1
                                end if
                                if(fx2<fm) then
                                        minimum = bx2
                                        fm = fx2
                                end if
                                if(fx2<fx1) then
                                        t = (bx1+bx2)/2.
                                        if(abs(bx1-t) < eps) return
                                        bx1=t
                                else
                                        t = (bx1+bx2)/2.
                                        if(abs(bx1-t) < eps)  return
                                        bx2=t
                                end if
                        end do
                end function


                real recursive function evaluate(self,x)
                        real :: x
                        class(funktion) :: self
                        !n=size(self%p_f_array)
                        evaluate = self%scal * self%p_eval(x)
                end function

                function evaluate_a(self,x)
                        real, dimension(:) :: x
                        real :: evaluate_a
                        class(funktion) :: self
                        !n=size(self%p_f_array)
                        evaluate_a = self%scal * self%p_eval_a(x)
                end function

                function init_dx()
                        type(funktion) :: init_dx
                        procedure(func), pointer :: p_id => id
                        init_dx%p_f =>  p_id !dx alone?!?
                        init_dx%p_div => dx_div
                        init_dx%p_mult => dx_mult
                end function        
                
                function init_x()
                        type(funktion) :: init_x
                        procedure(func), pointer :: p_id => id
                        init_x = create_funktion(p_id)
                end function
                function init_d()
                        type(funktion) :: init_d
                        procedure(func), pointer :: p_d => idd
                        init_d = create_funktion(p_d)
                end function
                real function idd(x)
                        real ::x
                        idd =1.0
                end function
                real function id(x)
                        real :: x
                        id =x
                end function
                subroutine init()
                        x = init_x()
                        d = init_d()
                        dx = init_dx()
                end subroutine
END MODULE
