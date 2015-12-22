!*****************************************************************************************
!>
!  Test program for [[ddeabm_class]].
!  Integrate a two-body orbit around the Earth.

    program test

    use ddeabm_module
    use root_module

    implicit none

    call zeroin_test()

    call ddeabm_test()

    end program test
!*****************************************************************************************
