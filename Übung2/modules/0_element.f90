module mElement
    implicit none
    type Element
        character(len=3) :: n
        real*8 :: pos(3)
    end type
end module
