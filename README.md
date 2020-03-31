# Experiments and performance evaluation of stochastic rounding algorithms
This repository countains the source code for reproducing the results in [Sec. 5, 1] and [Sec. 6, 1].

The performance benchmarks in `./performance` require the [GNU MPFR Library](https://mpfr.org/), which in turn depends on the [GNU Multiprecision Arithmetic (GMP) Library ](https://gmplib.org/). The codes can be compiled with `make all`, which generates the two executables `test_validation` and `test_performance`. The former performs a validation of the functions defined in [`stochrnddw.c`](./performance/stochrnddw.c), the latter produces the data in [Table 5.1, 1]. Compiling and running both tests is as easy as:
```console
cd performance
make all
make run_validation
make run_performance
```

The codes in `./numerical_experiments` generate the data used for the plots in [Fig. 6.1-6.5, 1]. The folder contains the following scripts:
   * [`test_summation_algorithms.m`](./numerical_experiments/test_summation_algorithms.m) [Sec. 6.1-6.2, 1].
   * [`test_summation_algorithms_double.m`](./numerical_experiments/test_summation_algorithms_double.m) [Sec. 6.2, 1].
   * [`ODE_tests.m`](./numerical_experiments/ODE_tests.m) [Sec. 6.3.1, 1].
   * [`unit_circle_ODE.m`](./numerical_experiments/unit_circle_ODE.m) [Sec. 6.3.2, 1].

The code for the experiment with double-precision arithmetic ([`test_summation_algorithms_double.m`](./numerical_experiments/test_summation_algorithms_double.m)) requires the [Stochastic Rounding Toolbox](https://github.com/mfasi/srtoolbox) be in the MATLAB path.

### References

 [1] M. Fasi and M. Mikaitis. [*Algorithms for stochastically rounded elementary arithmetic operations in IEEE 754 floating-point arithmetic*](http://eprints.maths.manchester.ac.uk/2758/), Technical Report 2020.9, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, Mar 2020.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENSE](./LICENSE)).
