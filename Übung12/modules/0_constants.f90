module constants
    character :: nl = char(10)
    real,parameter :: kb = 1.38064852e-23 ! J/K
    real,parameter :: e = 1.6021766208e-19 ! C
    real,parameter :: pi = 3.14159265359

    contains
    real function to(v, a,b)
       real :: v,a,b 
       to = v * a / b 
    end function
end module
