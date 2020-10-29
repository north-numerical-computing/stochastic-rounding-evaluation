#include "stochrnddw.h"

/* Compute stochastic rounding of a+b (new method). */
FPTYPE sr_add_rn_only(const FPTYPE a, const FPTYPE b) {
  // Compute floating-point approximation of sum and error.
  FPTYPE s, t;
  two_sum(&s, &t, a, b);

  // Check whether SR rounds up or down.
  FPTYPE round = sr_round(s, t);

  // Compute and return result.
  FPTYPE r = round + s;
  return r;
}

/* Compute stochastic rounding of a*b (new method). */
FPTYPE sr_mul_fma_rn_only(const FPTYPE a, const FPTYPE b) {
  // Compute floating-point approximation of product and error.
  FPTYPE s, t;
  two_prod_fma(&s, &t, a, b);

  // Check whether SR rounds up or down.
  FPTYPE round = sr_round(s, t);

  // Compute and return result.
  FPTYPE r = round + s;
  return r;
}

/* Compute stochastic rounding of a/b (new method). */
FPTYPE sr_div_rn_only(const FPTYPE a, const FPTYPE b) {
  // Compute floating-point quotient, remainder, and residual.
  FPTYPE s, t;
  s = a / b;
  t = fma(-s, b, a);
  t = t / b;

  // Check whether SR rounds up or down.
  FPTYPE round = sr_round(s, t);

  // Compute and return result.
  FPTYPE r = round + s;
  return r;
}

/* Compute stochastic rounding of sqrt(a) (new method). */
FPTYPE sr_sqrt_rn_only(const FPTYPE a) {
  FPTYPE s, t;
  s = sqrt(a);
  t = fma(-s, s, a);
  t = t/(2*s);

  // Check whether SR rounds up or down.
  FPTYPE round = sr_round(s, t);

  // Compute and return result.
  FPTYPE r = round + s;
  return r;
}
