#include <string.h>

#include "stochrndhp.h"
#include "stochrnddw.h"

/* Structure used for reference solution and probabilities. */
struct valinfo{
  double xd;       // RD(x) to double precision
  double rerrd;    // | x - xd | in double precision
  double probd;    // Probability of rounding towards xd
  double xu;       // RU(x) to double precision
  double rerru;    // | x - xu | in double precision
  double probu;    // Probability of rounding towards xu
  double epsilon;  // | xd - xu | in double precision
};

/* Populate valinfo structure from an MPFR value. */
void genvalinfo_mpfr (const mpfr_t x, struct valinfo* results, mpfr_t work) {
  // The variable work and results must be initialized outside this function!
  results->xd = mpfr_get_d(x, MPFR_RNDD);
  mpfr_set_d(work, results->xd, MPFR_RNDN);
  mpfr_sub(work, work, x, MPFR_RNDN);
  mpfr_abs(work, work, MPFR_RNDN);
  results->rerrd = mpfr_get_d(work, MPFR_RNDN);
  results->xu = mpfr_get_d(x, MPFR_RNDU);
  mpfr_set_d(work, results->xu, MPFR_RNDN);
  mpfr_sub(work, work, x, MPFR_RNDN);
  mpfr_abs(work, work, MPFR_RNDN);
  results->rerru = mpfr_get_d(work, MPFR_RNDN);
  results->epsilon = fabs(results->xu - results->xd);
  results->probd = results->rerru / results->epsilon;
  results->probu = results->rerrd / results->epsilon;
}

/* Populate valinfo structure from double precision roundings. */
double add_flop (const double a, const double b){return a+b;}
double sub_flop (const double a, const double b){return a-b;}
double mul_flop (const double a, const double b){return a*b;}
double div_flop (const double a, const double b){return a/b;}
void genvalinfo_op (double a, double b,
                    double (* op)(const double, const double),
                    int (* ref_op)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                    struct valinfo* results, mpfr_t work1, mpfr_t work2) {

  // Compute reference solution.
  mpfr_set_d(work1, a, MPFR_RNDN);
  mpfr_set_d(work2, b, MPFR_RNDN);
  ref_op(work1, work1, work2, MPFR_RNDN);

  // Rounding down.
  fesetround(FE_DOWNWARD);
  results->xd = op(a,b);
  mpfr_set_d(work2, results->xd, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerrd = mpfr_get_d(work2, MPFR_RNDN);

  // Rounding up.
  fesetround(FE_UPWARD);
  results->xu = op(a,b);
  mpfr_set_d(work2, results->xu, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerru = mpfr_get_d(work2, MPFR_RNDN);

  // Compute probabilities.
  results->epsilon = fabs(results->xu - results->xd);
  results->probd = results->rerru / results->epsilon;
  results->probu = results->rerrd / results->epsilon;
}

/* Generate pairs of random numbers and perform validation. */
void randtest_srop(size_t ntests, size_t nreps,
                   int (* mpfr_fun)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
                                    mpfr_rnd_t),
                   double (* test_fun)(double, double),
                   double (* sr_fun)(double, double)) {
  size_t i, j;
  double x, a, b;
  mpfr_t xh, ah, bh, work1, work2;
  mpfr_inits(xh, ah, bh, work1, work2, (mpfr_ptr) 0);
  struct valinfo results = {};
  size_t *counters = malloc(2*sizeof(size_t)); // counters[0]=#xd, counter[1]=#xu
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX*2.0-1.0;
    b = (double)rand()/RAND_MAX*2.0-1.0;
    mpfr_set_d(ah, a, MPFR_RNDN);
    mpfr_set_d(bh, b, MPFR_RNDN);
    mpfr_fun(xh, ah, bh, MPFR_RNDN); // reference solution in high precision
    genvalinfo_op(a, b, test_fun, mpfr_fun, &results, work1, work2);
    memset(counters, 0, 2*sizeof(size_t));
    for (j=0; j<nreps; j++) {
      x = sr_fun(a, b);
      if (x == results.xd)
        counters[0]++;
      else if (x == results.xu)
        counters[1]++;
      else {
        printf("\n\n[ERROR] Not rounding to either nearest fp value.\n");
        exit(10);
      }
    }
    if(fabs(counters[0]/(double)nreps - results.probd) > 1e-2) {
      printf("\n\n[ERROR] The test failed for a = %.16e, b = %.16e.\n", a, b);
      exit(1);
    }
  }
}


int main (int argc, char *argv[]) {
  // Initialize MPFR and PRNGs.
  size_t i, j;
  const int p = 113; // Precision bits in quad precision.
  const unsigned long seed = 1;
  static gmp_randstate_t randstate;
  gmp_randinit_default(randstate);
  gmp_randseed_ui(randstate, seed);
  mpfr_set_default_prec(p);
  srand(seed);

  // Validation initializations.
  mpfr_t xh, xdmp, work;
  mpfr_inits2(p, xh, xdmp, work, (mpfr_ptr) 0);
  struct valinfo results = {};
  double x;
  size_t *counters = malloc(2*sizeof(size_t)); // counters[0]=#xd, counter[1]=#xu
  size_t ntests = 1000, nreps = 50000;

  // Validate stochastic rounding function sr_mpfr().
  printf("[INFO] Now testing sr_mpfr()... ");
  for (i=0; i<ntests; i++) {
    mpfr_urandom(xh, randstate, MPFR_RNDN);
    genvalinfo_mpfr(xh, &results, work);
    memset(counters, 0, 2*sizeof(size_t));
    for (j=0; j<nreps; j++) {
      x = sr_mpfr(xh, xdmp, work);
      if (x == results.xd)
        counters[0]++;
      else if (x == results.xu)
        counters[1]++;
      else {
        printf("[ERROR] Not rounding to either nearest fp value.");
        exit(10);
      }
    }
    if(fabs(counters[0]/(double)nreps - results.probd) > 1e-1) {
      printf("[ERROR] The test sr_mpfr() failed for xh = ");
      mpfr_out_str(stdout, 10, 0, xh, MPFR_RNDD);
      printf(".\n");
      exit(1);
    }
  }
  printf("all tests passed.\n");

  // Validate sr_add().
  printf("[INFO] Now testing sr_add()... ");
  randtest_srop(ntests, nreps, mpfr_add, add_flop, sr_add);
  printf("all tests passed.\n");

  // Validate fast_sr_add().
  printf("[INFO] Now testing fast_sr_add()... ");
  randtest_srop(ntests, nreps, mpfr_add, add_flop, fast_sr_add);
  printf("all tests passed.\n");

  // Validate sr_mul_fma().
  printf("[INFO] Now testing sr_mul_fma()... ");
  randtest_srop(ntests, nreps, mpfr_mul, mul_flop, sr_mul_fma);
  printf("all tests passed.\n");

  // Validate sr_div().
  printf("[INFO] Now testing sr_div()... ");
  randtest_srop(ntests, nreps, mpfr_div, div_flop, sr_div);
  printf("all tests passed.\n");

  // Terminate.
  free(counters);
  mpfr_clears(xh, xdmp, work, (mpfr_ptr) 0);
  mpfr_free_cache ();

  return 0;
}
