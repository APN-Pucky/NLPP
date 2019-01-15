module analysis
    use constants, only: pi
    use parameters, only: N,dt,max_iter,eps,h
    use funktions


    interface nst_newton
        module procedure nst_newton_f,nst_newton_fdf,nst_bisection
    end interface
 
    interface dft
        module procedure dft_f,dft_y
    end interface
    interface idft
        module procedure idft_y!,idft_f
    end interface

    contains
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

                function dft_y(dy,t_start)
                        complex :: dy(N)
                        complex :: dft_y(N+1)
                        complex :: i = (0.,1.)
                        integer :: nn,k
                        real :: t_start
                        dft_y = (/(sum((/(dy(k+1)*exp(-i*2*pi/N/dt*(t_start+dt*k)*nn),k=0,N-1)/)),nn=0,N)/)
                end function

                function idft_y(dy,t_start)
                        complex :: dy(N+1)
                        complex :: idft_y(N)
                        complex :: i = (0.,1.)
                        integer :: nn,k
                        real :: t_start
                        idft_y = (/(sum((/(dy(nn+1)*exp(i*2*pi/N/dt*(t_start+dt*k)*nn),nn=0,N-1)/)),k=0,N)/)/N
                end function

                function dft_f(f,t_start) !set dt and N before
                        complex :: dft_f(N+1)
                        type(funktion) :: f
                        integer :: i
                        real :: t_start
                        dft_f = dft_y((/(cmplx(f%get(t_start+dt*i),0.),i=0,N-1)/),t_start)
                end function

end module
