module mMath
    contains
    real function average(array)
        real*8 :: array(:)
        average = SUM(array)/size(array)
    end function
    real function variance(array)
        real*8 :: array(:)
        real*8 :: tmp(size(array))
        variance = -average(array)**2
        tmp = array**2
        variance = variance + average(tmp)
    end function
    subroutine math_test
        real*8 :: arr(4)
        arr(1) = 1
        arr(2) = 2
        arr(3) = 3
        arr(4) = 4
        !arr(4) = variance(arr)
        write(*,*) average(arr), variance(arr)
    end subroutine
end module
