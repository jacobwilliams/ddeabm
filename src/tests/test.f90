!*****************************************************************************************
!>
!  Test program for ddeabm.

    program test

    use ddeabm_module, only: ddeabm_test
    use root_module,   only: zeroin_test

    implicit none

    call zeroin_test()
    call ddeabm_test()

    end program test
!*****************************************************************************************
