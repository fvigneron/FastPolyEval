# Attribution notice

## Included libraries

- This software includes the [MPFR](https://www.mpfr.org) library which has his own separate license
(see [GNU LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.html)) to handle arbitrary high-precision arithmetic.
- In turn, MPFR includes the [GMP](https://gmplib.org) library which can either be distributed
under [GNU LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.html) or [GNU GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

## Alternative libraries

It is possible to compile a version of the code that uses hardware arithmetic at runtime, up to FP80 numbers.
For now, the MPFR library is included by default, even if you inted to use only FP80 arithmetic.

We welcome the idea of using alternative libraries for high-precision arithmetic, in particular if your
favorite library works on GPU or mixed CPU/GPU environment. If you are interested, please contact me
(<francois.vigneron@univ-reims.fr>) so we can discuss evolving our software accordingly :wink:.
