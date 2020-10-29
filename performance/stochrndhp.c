#include "stochrndhp.h"
#include "stochrnddw.h"

/* Compute stochastic rounding of x using MPFR. */
double sr_mpfr(mpfr_t x, mpfr_t w1, mpfr_t w2) {
  // Compute sign and truncation of absolute value in double precision.
  double signx = mpfr_signbit(x)==1? -1: 1; // sign(x)

  double xd = mpfr_get_d(x, MPFR_RNDZ);
  mpfr_set_d(w1, xd, MPFR_RNDZ);

  // Compute absolute value of residual.
  mpfr_sub(w1, x, w1, MPFR_RNDZ);
  mpfr_abs(w1, w1, MPFR_RNDZ);
  xd = fabs(xd);

  // Rescale residual between 0 and 1.
  int exponent = get_exponent(xd);
  mpfr_mul_2exp(
    w1, w1, 52-exponent, MPFR_RNDN); // hardcoded for double precision

  // Perform stochastic rounding.
  double rnd = (double)rand() / RAND_MAX;
  mpfr_set_d(w2, rnd, MPFR_RNDN);

  return signx * (mpfr_lessequal_p(w2,w1)? nextafter(xd,INFINITY) : xd);
}


/* Compute stochastic rounding of op(a,b) using MPFR. */
double sr_op_mpfr_noalloc(const double a, const double b,
                          int (* mpfr_op)(mpfr_ptr,
                                          mpfr_srcptr, mpfr_srcptr,
                                          mpfr_rnd_t),
                          mpfr_t w1, mpfr_t w2, mpfr_t w3) {
  mpfr_set_d(w1, a, MPFR_RNDN);
  mpfr_set_d(w2, b, MPFR_RNDN);
  mpfr_op(w1, w1, w2, MPFR_RNDN);
  return sr_mpfr(w1, w2, w3);
}

double sr_op_mpfr_noalloc_single_arg(const double a,
                          int (* mpfr_op)(mpfr_ptr,
                                          mpfr_srcptr,
                                          mpfr_rnd_t),
                          mpfr_t w1, mpfr_t w2, mpfr_t w3) {
  mpfr_set_d(w1, a, MPFR_RNDN);
  mpfr_op(w1, w1, MPFR_RNDN);
  return sr_mpfr(w1, w2, w3);
}

/* Compute stochastic rounding of a+b using MPFR. */
double sr_add_mpfr(const double a, const double b) {
  // Initialize parameters.
  mpfr_t w1, w2, w3;
  mpfr_inits(w1, w2, w3, (mpfr_ptr) 0);
  double z = sr_op_mpfr_noalloc(a, b, mpfr_add, w1, w2, w3);
  mpfr_clears(w1, w2, w3, (mpfr_ptr) 0);
  return z;
}

/* Compute stochastic rounding of a-b using MPFR. */
double sr_sub_mpfr(const double a, const double b) {
  // Initialize parameters.
  mpfr_t w1, w2, w3;
  mpfr_inits(w1, w2, w3, (mpfr_ptr) 0);
  double z = sr_op_mpfr_noalloc(a, b, mpfr_sub, w1, w2, w3);
  mpfr_clears(w1, w2, w3, (mpfr_ptr) 0);
  return z;
}

/* Compute stochastic rounding of a*b using MPFR. */
double sr_mul_mpfr(const double a, const double b) {
  mpfr_t w1, w2, w3;
  mpfr_inits(w1, w2, w3, (mpfr_ptr) 0);
  double z = sr_op_mpfr_noalloc(a, b, mpfr_mul, w1, w2, w3);
  mpfr_clears(w1, w2, w3, (mpfr_ptr) 0);
  return z;
}

/* Compute stochastic rounding of a/b using MPFR. */
double sr_div_mpfr(const double a, const double b) {
  mpfr_t w1, w2, w3;
  mpfr_inits(w1, w2, w3, (mpfr_ptr) 0);
  double z = sr_op_mpfr_noalloc(a, b, mpfr_div, w1, w2, w3);
  mpfr_clears(w1, w2, w3, (mpfr_ptr) 0);
  return z;
}

/* Compute stochastic rounding of sqrt(a) using MPFR. */
double sr_sqrt_mpfr(const double a) {
  mpfr_t w1, w2, w3;
  mpfr_inits(w1, w2, w3, (mpfr_ptr) 0);
  double z = sr_op_mpfr_noalloc_single_arg(a, mpfr_sqrt, w1, w2, w3);
  mpfr_clears(w1, w2, w3, (mpfr_ptr) 0);
  return z;
}
