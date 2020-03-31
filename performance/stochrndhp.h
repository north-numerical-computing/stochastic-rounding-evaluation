#ifndef STOCHRNDHP_H_
#define STOCHRNDHP_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gmp.h>
#include <mpfr.h>

double sr_mpfr(mpfr_t, mpfr_t, mpfr_t);
double sr_op_mpfr_noalloc(const double, const double,
                          int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                          mpfr_t, mpfr_t, mpfr_t);

double sr_add_mpfr(const double, const double);
double sr_sub_mpfr(const double, const double);
double sr_mul_mpfr(const double, const double);
double sr_div_mpfr(const double, const double);


#endif // STOCHRNDHP_H_
