#ifndef STOCHRNDDW_H_
#define STOCHRNDDW_H_

#include <assert.h>
#include <stdlib.h>
#include <fenv.h>
#include <math.h>

#define SIGN(x) ((x>0) - (x<0))

void two_sum(double*, double*, const double, const double);
void fast_two_sum(double*, double*, const double, const double);
void two_prod_fma(double*, double*, const double, const double);

double sr_add(const double, const double);
double fast_sr_add(const double, const double);
double sr_mul_fma(const double, const double);
double sr_div(const double, const double);

#endif // STOCHRNDDW_H_
