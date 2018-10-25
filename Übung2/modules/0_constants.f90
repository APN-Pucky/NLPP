module mConstants
    character :: nl = char(10)
    real*8,parameter :: kb = 1.38064852e-23 ! J/K
    real*8,parameter :: e = 1.6021766208e-19 ! C

    contains
    real*8 function to(v, a,b)
       real*8 :: v,a,b 
       to = v * a / b 
    end function
end module
