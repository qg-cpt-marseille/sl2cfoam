## Numerical methods for the EPRL spin foam transition amplitudes and Lorentzian recoupling theory

The intricated combinatorial structure and the non-compactness of the group have always made the computation of SL(2,C) EPRL spin foam transition amplitudes a very hard and resource demanding task. With sl2cfoam we provide a C-coded library for the evaluation of the Lorentzian EPRL vertex amplitude.
We provide a tool to compute the Lorentzian EPRL vertex amplitude in the intertwiner basis and some utilities to evaluate SU(2) invariants, booster functions and SL(2,C) Clebsch-Gordan coefficients.

### Compiling and using

Linux/Unix system supported. Few libraries are needed :

* Wigxjpf (http://fy.chalmers.se/subatom/wigxjpf/) to compute 3j,6j and 9j symbols. This is included in the `ext/` directory and must be built before compiling sl2cfoam.

* GMP (https://gmplib.org/) to compute d-small Wigner matrices for coherent states and for booster function.

* MPFR (http://www.mpfr.org). As GMP, high precision library for floating points. Used in coherent states and boosters computations.

* MPC (http://www.multiprecision.org/mpc/). As GMP and MPFR, high precision library for complex numbers. Used for booster functions.

* GSL (https://www.gnu.org/software/gsl/doc/html/intro.html) to compute complex gamma functions in booster functions.

* OpenMp (http://www.openmp.org/) to parallelize the computation.

To build the sl2cfoam library you can type:

* `make` creates the library and compiles the test programs.
* `make test` compiles only test programs.
* `./testall` launches test programs.

Then, to link the library to your program you can use the flag `-lsl2cfoam`. Also, you have to include the main header `inc/sl2cfoam.h`.
In `test/` there are several test files that you can look at to learn how to use the library. 

### Licenses

sl2cfoam is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

sl2cfoam is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
