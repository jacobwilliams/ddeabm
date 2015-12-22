!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!
!  Numeric kind definitions.

    module kind_module

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    integer,parameter,public :: wp = real64  !! Using "double precision" real kinds

    end module kind_module
!*****************************************************************************************
