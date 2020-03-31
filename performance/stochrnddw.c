#include "stochrnddw.h"

/* Compute s and t such that s+t = a+b. */
void two_sum(double* s, double* t,
             const double a, const double b) {
  *s = a + b;
  double ap = *s - b;
  double bp = *s - ap;
  ap = a - ap;
  bp = b - bp;
  *t = ap + bp;
}

/* If |a| > |b|, compute s and t such that s+t = a+b. */
void fast_two_sum(double* s, double* t,
                  const double a, const double b) {
  assert(fabs(a) > fabs(b));
  *s = a + b;
  *t = (*s-a) - b;
}

/* Compute s and t such that s+t = a*b. */
void two_prod_fma(double* s, double* t,
                  const double a, const double b) {
  *s = a * b;
  *t = fma(a, b, -*s);
}

/* Compute stochastic rounding of a+b (new method). */
double sr_add(const double a, const double b) {
  int prevround = fegetround();
  fesetround(FE_TONEAREST);

  // Compute floating-point approximation of sum and error.
  double s, t;
  two_sum(&s, &t, a, b);

  // Compute exponent with truncation.
  int exponent;
  fesetround(FE_TOWARDZERO);
  frexp(a+b,&exponent); // Exponent of fraction / 2 (i.e. in [0.5, 1)).
  /* assert(exponent-1 == floor(log2(fabs(a+b)))); */

  // Compute and renormalize random number.
  double p = ldexp(SIGN(t) * (rand()/(double)RAND_MAX), exponent-53);

  // Compute and return result.
  if (t >= 0)
    fesetround(FE_DOWNWARD);
  else
    fesetround(FE_UPWARD);
  double r = (t+p) + s;
  fesetround(prevround);
  return r;
}

/* Compute approximate stochastic rounding of a+b (new method) */
double fast_sr_add(const double a, const double b) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point approximation of sum and error.
  double s, t;
  two_sum(&s, &t, a, b);

  // Compute exponent with truncation.
  int exponent;
  frexp(s, &exponent);

  // Compute and renormalize random number.
  double p = ldexp(SIGN(t) * (rand()/(double)RAND_MAX), exponent-53);

  // Compute and return result.
  double r = (t+p) + s;
  fesetround(prevround);
  return r;
}

/* Compute stochastic rounding of a*b (new method). */
double sr_mul_fma(const double a, const double b) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point approximation of product and error.
  double s, t;
  two_prod_fma(&s, &t, a, b);

  // Compute exponent with truncation.
  int exponent;
  frexp(s, &exponent);

  // Compute and renormalize random number.
  double p = ldexp(SIGN(t) * (rand()/(double)RAND_MAX), exponent-53);

  // Compute and return result.
  double r = (t+p) + s;
  fesetround(prevround);
  return r;
}

/* Compute stochastic rounding of a/b (new method). */
double sr_div(const double a, const double b) {
  int prevround = fegetround();
  fesetround(FE_TOWARDZERO);

  // Compute floating-point quotient, remainder, and residual.
  double s, t;
  s = a / b;
  t = fma(-s, b, a);
  t = t / b;

  // Compute exponent with truncation.
  int exponent;

  frexp(s, &exponent);

  // Compute and renormalize random number.
  double p = ldexp(SIGN(t) * (rand()/(double)RAND_MAX), exponent-53);

  // Compute and return result.
  double r = (t+p) + s;
  fesetround(prevround);
  return r;
}
