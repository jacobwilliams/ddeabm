project: ddeabm
project_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/ddeabm
summary: Modern Fortran Implementation of the DDEABM Adams-Bashforth-Moulton ODE Solver
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
exclude_dir: ./tests
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran
            roots-fortran:https://github.com/jacobwilliams/roots-fortran
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}
