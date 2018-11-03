program main
    use io
    use random
    implicit none
    call io_selectLoop((/&
        Select("Set RNG Mode", set_rng_mode),&
        Select("Set Fortran Mode", set_fortran_mode),&
        Select("Next RNG Number", random_print_next),&
        Select("Save n Numbers", random_save_next)&
    /))

end program
