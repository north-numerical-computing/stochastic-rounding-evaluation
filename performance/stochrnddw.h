#ifndef STOCHRNDDW_H_
#define STOCHRNDDW_H_

#include <assert.h>
#include <stdlib.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>

#define FPTYPE    double
#define INTTYPE uint64_t
#define SINTTYPE int64_t
#define DEFPREC       53
#define EMAX        1023
#define NBITS         64
#define NLEADBITS     12
#define NRNDBITS      30
#define FULLMASK  0xFFFFFFFFFFFFFFFF
#define ABSMASK   0x7FFFFFFFFFFFFFFF
#define EXPMASK   0x7FF0000000000000
#define FRACMASK  0x000FFFFFFFFFFFFF
#define SIGNMASK  0x8000000000000000

#define INTCONST(x)(INTTYPE)(x)

typedef union {
  FPTYPE fpval;
  INTTYPE intval;
} STRUCTNAME;

#define INTOF(x)(((STRUCTNAME *)(x))->intval)
#define FPOF(x)(((STRUCTNAME)(x)).fpval)


#define SIGN(x) ((x>0) - (x<0))

#define FASTSIGN(x)(SIGNMASK & INTOF(x))
#define INTABS(x)(ABSMASK & INTOF(x))
#define ABS(x)(FPOF(ABSMASK & INTOF(x)))

/************************
 * AUGMENTED OPERATIONS *
 ************************/

/* Compute s and t such that s+t = a+b. */
static inline void two_sum(FPTYPE* s, FPTYPE* t,
                           const FPTYPE a, const FPTYPE b) {
  *s = a + b;
  FPTYPE ap = *s - b;
  FPTYPE bp = *s - ap;
  ap = a - ap;
  bp = b - bp;
  *t = ap + bp;
}

/* If |a| > |b|, compute s and t such that s+t = a+b. */
static inline void fast_two_sum(FPTYPE* s, FPTYPE* t,
                                const FPTYPE a, const FPTYPE b) {
  assert(fabs(a) > fabs(b));
  *s = a + b;
  *t = (*s-a) - b;
}

/* Compute s and t such that s+t = a*b. */
static inline void two_prod_fma(FPTYPE* s, FPTYPE* t,
                                const FPTYPE a, const FPTYPE b) {
  *s = a * b;
  *t = fma(a, b, -*s);
}

/*********************
 * UTILITY FUNCTIONS *
 *********************/

/* Compute IEEE floating-point exponent of x. */
static inline int get_exponent(FPTYPE x) {
  return (((INTOF(&x) & EXPMASK) >> 52) - 1023);
}

/* Compute absolute value of ulp(x). */
static inline INTTYPE ulp_abs(FPTYPE x) {
  SINTTYPE tmp = get_exponent(x);
  if (tmp <= (-1023+52)) // ulp(x) is subnormal
    return INTCONST(1);
  else                   // ulp(x) i normal
    return INTCONST(tmp + (+1023-52)) << 52;
}

/* A helper function for SR algorithms without the change of rounding mode. */
static inline FPTYPE sr_round(const FPTYPE s, const FPTYPE t) {
  FPTYPE ulp;
  if (FASTSIGN(&t) != FASTSIGN(&s))
    //    ulp = FPOF(FASTSIGN(&t) | ulp_abs(s*(1 - ldexp(1, -53))));
    ulp = SIGN(t) * ldexp(1, get_exponent(s*(1-ldexp(1, -53)))-52);
  else
    //    ulp = FPOF(FASTSIGN(&t) | ulp_abs(s));
    ulp = SIGN(t) * ldexp(1, get_exponent(s)-52);

  FPTYPE p = (rand()/(FPTYPE)RAND_MAX) * ulp;
  FPTYPE round = 0;
  if (fabs(t+p) >= fabs(ulp))
    round = ulp;
  return round;
}

/****************************************************
 * ALGORITHMS WITH EXPLICIT CHANGE OF ROUNDING MODE *
 ****************************************************/

FPTYPE sr_add(const FPTYPE, const FPTYPE);
FPTYPE sr_mul_fma(const FPTYPE, const FPTYPE);
FPTYPE sr_div(const FPTYPE, const FPTYPE);
FPTYPE sr_sqrt(const FPTYPE);

/*****************************************************
 * ALGORITHMS WITH SIMULATED CHANGE OF ROUNDING MODE *
 *****************************************************/

FPTYPE sr_add_rn_only(const FPTYPE a, const FPTYPE b);
FPTYPE sr_mul_fma_rn_only(const FPTYPE, const FPTYPE);
FPTYPE sr_div_rn_only(const FPTYPE, const FPTYPE);
FPTYPE sr_sqrt_rn_only(const FPTYPE);

#endif // STOCHRNDDW_H_
