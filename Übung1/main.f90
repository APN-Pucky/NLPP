program main
    use pi
    use sampling
    implicit none
    write(*,100) "Pi"
    call pi_main()
    write(*,100) "Monte-Carlo Sampling"
    call sampling_main()

    !formats!
    100 FORMAT(50("="),T8,A)
    101 FORMAT(50("-"),T5,A)
end program



