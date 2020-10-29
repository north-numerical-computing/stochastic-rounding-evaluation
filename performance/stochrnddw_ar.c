#include "stochrnddw.h"

/* Compute stochastic rounding of a+b. */
FPTYPE sr_add(const FPTYPE a, const FPTYPE b) {
  int prevround = fegetround();
  fesetround(FE_TONEAREST);

  // Compute floating-point approximation of sum and error.
  FPTYPE s, t;
  two_sum(&s, &t, a, b);

  // Compute exponent with truncation.
  int exponent;
  fesetround(FE_TOWARDZERO);
  exponent = get_exponent(a+b); // Exponent of fraction / 2 (i.e. in [0.5, 1)).

  // Compute and renormalize a random number.
  FPTYPE p = ldexp(SIGN(t) * (rand()/(FPTYPE)RAND_MAX), exponent-52);

  // Compute and return result.
  if (t >= 0)
    fesetround(FE_DOWNWARD);
  else
    fesetround(FE_UPWARD);

  FPTYPE r = (t + p) + s;

  fesetround(prevround);
  return r;
}

/* Compute stochastic rounding of a*b. */
FPTYPE sr_mul_fma(const FPTYPE a, const FPTYPE b) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point approximation of product and error.
  FPTYPE s, t;
  two_prod_fma(&s, &t, a, b);

  // Compute exponent with truncation.
  int exponent = get_exponent(s);

  // Compute and renormalize random number.
  FPTYPE p = ldexp(SIGN(t) * (rand()/(FPTYPE)RAND_MAX), exponent-52);

  // Compute and return result.
  FPTYPE r = (t + p) + s;

  fesetround(prevround);
  return r;
}

/* Compute stochastic rounding of a/b. */
FPTYPE sr_div(const FPTYPE a, const FPTYPE b) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point quotient, remainder, and residual.
  FPTYPE s, t;
  s = a / b;
  t = fma(-s, b, a);
  t = t / b;

  // Compute exponent with truncation.
  int exponent = get_exponent(s);

  // Compute and renormalize random number.
  FPTYPE p = ldexp(SIGN(t) * (rand()/(FPTYPE)RAND_MAX), exponent-52);

  // Compute and return result.
  FPTYPE r = (t + p) + s;

  fesetround(prevround);
  return r;
}

/* Compute stochastic rounding of sqrt(a). */
FPTYPE sr_sqrt(const FPTYPE a) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point quotient, remainder, and residual.
  FPTYPE s, t;
  s = sqrt(a);
  t = fma(-s, s, a);
  t = t/(2*s);

  // Compute exponent with truncation.
  int exponent = get_exponent(s);

  // Compute and renormalize random number.
  FPTYPE p = ldexp(SIGN(t) * (rand()/(FPTYPE)RAND_MAX), exponent-52);

  // Compute and return result.
  FPTYPE r = (t + p) + s;

  fesetround(prevround);
  return r;
}
