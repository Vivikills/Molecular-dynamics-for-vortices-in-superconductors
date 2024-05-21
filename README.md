Files represent a molecular dynamics project for vortices in a superconductor which was written in Fortran. Short description of the subroutines is below:
- rand1.f90 generates random coordinates of vortices in the boundary region when they first enter the sample
- processes1.f90 is responsible for the processes occurring in the system, such as the birth and destruction of vortexes, calculation of energies, etc.
- LJ_DE10.f90 is responsible for launching the entire algorithm, that is, for solving a system of differential equations using the Runge-Kutta method
