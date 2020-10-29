#include <string.h>
#include <float.h>

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
double sqrt_flop (const double a){return sqrt(a);}

void genvalinfo_op (double a, double b,
                    double (* op)(const double, const double),
                    int (* ref_op)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t),
                    struct valinfo* results, mpfr_t work1, mpfr_t work2) {

  // Compute a reference solution.
  mpfr_set_d(work1, a, MPFR_RNDN);
  mpfr_set_d(work2, b, MPFR_RNDN);
  ref_op(work1, work1, work2, MPFR_RNDN);

  // Rounding down.
  int prevround = fegetround();
  fesetround(FE_DOWNWARD);
  results->xd = op(a, b);
  mpfr_set_d(work2, results->xd, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerrd = mpfr_get_d(work2, MPFR_RNDN);

  // Rounding up.
  fesetround(FE_UPWARD);
  results->xu = op(a, b);
  mpfr_set_d(work2, results->xu, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerru = mpfr_get_d(work2, MPFR_RNDN);
  fesetround(prevround);

  // Compute probabilities.
  results->epsilon = fabs(results->xu - results->xd);
  results->probd = results->rerru / results->epsilon;
  results->probu = results->rerrd / results->epsilon;
}

void genvalinfo_op_single_arg (double a,
                               double (* op)(const double),
                               int (* ref_op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                               struct valinfo* results, mpfr_t work1, mpfr_t work2) {

  // Compute a reference solution.
  mpfr_set_d(work1, a, MPFR_RNDN);
  ref_op(work1, work1, MPFR_RNDN);

  // Rounding down.
  int prevround = fegetround();
  fesetround(FE_DOWNWARD);
  results->xd = op(a);
  mpfr_set_d(work2, results->xd, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerrd = mpfr_get_d(work2, MPFR_RNDN);

  // Rounding up.
  fesetround(FE_UPWARD);
  results->xu = op(a);
  mpfr_set_d(work2, results->xu, MPFR_RNDN);
  mpfr_sub(work2, work1, work2, MPFR_RNDN);
  mpfr_abs(work2, work2, MPFR_RNDN);
  results->rerru = mpfr_get_d(work2, MPFR_RNDN);
  fesetround(prevround);

  // Compute probabilities.
  results->epsilon = fabs(results->xu - results->xd);
  results->probd = results->rerru / results->epsilon;
  results->probu = results->rerrd / results->epsilon;
}

void perform_a_test(double a, double b, double expected1, double expected2,
                    size_t nreps,
                    double prob, // Prob. of obraining the value expected1
                    double (* sr_fun)(double, double)) {

  double x;
  size_t *counters = malloc(2*sizeof(size_t));
  counters[0] = 0; counters[1] = 0;
  memset(counters, 0, 2*sizeof(size_t));
  for (size_t j=0; j<nreps; j++) {
    x = sr_fun(a, b);
    if (x == expected1)
      counters[0]++;
    else if (x == expected2)
      counters[1]++;
    else {
      printf("\n\n[ERROR] Not rounding to either nearest fp values.\n");
      exit(10);
    }
  }
  if(fabs(counters[0]/(double)nreps - prob) > 1e-2) {
    printf("\n\n[ERROR] The test failed for a = %.16e, b = %.16e.\n", a, b);
    printf("[ERROR] Expected %.4f but obtained %.4f.\n",
           prob,
           counters[0]/(double)nreps);
    exit(1);
  }
}

void special_cases_srmul(size_t nreps,
                         double (* sr_fun)(double, double)) {

  perform_a_test(DBL_TRUE_MIN, 0.5, 0, DBL_TRUE_MIN, nreps, 1.0, sr_fun);
  perform_a_test(DBL_MAX, DBL_MAX, INFINITY, DBL_MAX, nreps, 0.0, sr_fun);
}

void special_cases_srmul_rn_only(size_t nreps,
                                 double (* sr_fun)(double, double)) {

  perform_a_test(DBL_TRUE_MIN, 0.5, 0, DBL_TRUE_MIN, nreps, 1.0, sr_fun);
  perform_a_test(DBL_MAX, DBL_MAX, INFINITY, DBL_MAX, nreps, 1.0, sr_fun);
}

/* Test some special cases */
void special_tests_sradd_rn_only(size_t nreps,
                                 double (* sr_fun)(double, double)) {

  // Test that we round with a correct prob. when x lands between
  // the maximum value and the threshold for rounding to inf.
  perform_a_test(DBL_MAX,
                 ldexp(0.25, 971), // Quarter of an ulp of DBL_MAX
                 DBL_MAX, INFINITY, nreps, 0.75, sr_fun);

  // Test that we return infinity in 100% cases when result crosses
  // the threshold.
  perform_a_test(DBL_MAX, DBL_MAX, DBL_MAX, INFINITY, nreps, 0.0, sr_fun);

  // Test a = -b, in which 0 should be always returned.
  perform_a_test(1, -1, 0, DBL_TRUE_MIN, nreps, 1.0, sr_fun);
}

/* Test some special cases */
void special_tests_sradd(size_t nreps,
                         double (* sr_fun)(double, double)) {

  perform_a_test(DBL_MAX, ldexp(0.25, 971), // Quarter of an ulp of DBL_MAX
                 DBL_MAX, INFINITY, nreps, 1.0, sr_fun);

  // Test a = -b, in which 0 should be always returned.
  perform_a_test(1, -1, 0, DBL_TRUE_MIN, nreps, 1.0, sr_fun);

  // Test that we return infinity in 100% cases when result crosses
  // the threshold.
  size_t *counters = malloc(2*sizeof(size_t));
  double x;
  double a = DBL_MAX;
  double b = DBL_MAX;
  counters[0] = 0; counters[1] = 0;
  for (size_t j=0; j<nreps; j++) {
    x = sr_fun(a, b);
    if (isnan(x))
      counters[0]++;
    else if (x == INFINITY)
      counters[1]++;
    else {
      printf("\n\n[ERROR] Not rounding to either nearest fp values. %f \n", x);
      exit(10);
    }
  }
  if(fabs(counters[0]/(double)nreps - 1) > 1e-2) {
    printf("\n\n[ERROR] The test failed for a = %.16e, b = %.16e.\n", a, b);
    printf("[ERROR] Expected %.4f but obtained %.4f.\n",
           1.0,
           counters[0]/(double)nreps);
    exit(1);
  }
}

/* Generate pairs of random numbers and perform validation. */
void randtest_srop(size_t ntests, size_t nreps,
                   int (* mpfr_fun)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
                                    mpfr_rnd_t),
                   double (* test_fun)(double, double),
                   double (* sr_fun)(double, double)) {

  size_t i;
  double a, b;
  mpfr_t xh, ah, bh, work1, work2;
  mpfr_inits(xh, ah, bh, work1, work2, (mpfr_ptr) 0);
  struct valinfo results = {};
  size_t *counters = malloc(2*sizeof(size_t));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX*2.0-1.0;
    b = (double)rand()/RAND_MAX*2.0-1.0;
    mpfr_set_d(ah, a, MPFR_RNDN);
    mpfr_set_d(bh, b, MPFR_RNDN);
    mpfr_fun(xh, ah, bh, MPFR_RNDN); // reference solution in high precision
    genvalinfo_op(a, b, test_fun, mpfr_fun, &results, work1, work2);
    memset(counters, 0, 2*sizeof(size_t));
    perform_a_test(a, b, results.xd, results.xu, nreps, results.probd, sr_fun);
  }
}

/* Generate pairs of random numbers and perform validation. */
void randtest_srop_single_arg(size_t ntests, size_t nreps,
                              int (* mpfr_fun)(mpfr_ptr, mpfr_srcptr,
                                               mpfr_rnd_t),
                              double (* test_fun)(double),
                              double (* sr_fun)(double)) {

  size_t i, j;
  double x, a;
  mpfr_t xh, ah, work1, work2;
  mpfr_inits(xh, ah, ah, work1, work2, (mpfr_ptr) 0);
  struct valinfo results = {};
  size_t *counters = malloc(2*sizeof(size_t));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX;
    mpfr_set_d(ah, a, MPFR_RNDN);
    mpfr_fun(xh, ah, MPFR_RNDN); // reference solution in high precision
    genvalinfo_op_single_arg(a, test_fun, mpfr_fun, &results, work1, work2);
    memset(counters, 0, 2*sizeof(size_t));
    for (j=0; j<nreps; j++) {
      x = sr_fun(a);
      if (x == results.xd)
        counters[0]++;
      else if (x == results.xu)
        counters[1]++;
      else {
        printf("\n\n[ERROR] Not rounding to either nearest fp values.");
        exit(10);
      }
    }
    if(fabs(counters[0]/(double)nreps - results.probd) > 1e-2) {
      printf("\n\n[ERROR] The test failed for a = %.16e. %e %e \n",
             a, counters[0]/(double)nreps, results.probd);
      exit(1);
    }
  }
}


int main () {

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
  size_t *counters = malloc(2*sizeof(size_t));
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
        printf("[ERROR] Not rounding to either nearest fp values.");
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


  /* Test SR operation with random inputs. */

  // Validate sr_add().
  printf("[INFO] Now testing sr_add()... ");
  randtest_srop(ntests, nreps, mpfr_add, add_flop, sr_add);
  printf("all tests passed.\n");

  // Validate sr_add_rn_only().
  printf("[INFO] Now testing sr_add_rn_only()... ");
  randtest_srop(ntests, nreps, mpfr_add, add_flop, sr_add_rn_only);
  printf("all tests passed.\n");

  // Validate sr_mul_fma().
  printf("[INFO] Now testing sr_mul_fma()... ");
  randtest_srop(ntests, nreps, mpfr_mul, mul_flop, sr_mul_fma);
  printf("all tests passed.\n");

  // Validate sr_mul_fma_rn_only().
  printf("[INFO] Now testing sr_mul_fma_rn_only()... ");
  randtest_srop(ntests, nreps, mpfr_mul, mul_flop, sr_mul_fma_rn_only);
  printf("all tests passed.\n");

  // Validate sr_div().
  printf("[INFO] Now testing sr_div()... ");
  randtest_srop(ntests, nreps, mpfr_div, div_flop, sr_div);
  printf("all tests passed.\n");

  // Validate sr_div_rn_only().
  printf("[INFO] Now testing sr_div_rn_only()... ");
  randtest_srop(ntests, nreps, mpfr_div, div_flop, sr_div_rn_only);
  printf("all tests passed.\n");

  // Validate sr_sqrt().
  printf("[INFO] Now testing sr_sqrt()... ");
  randtest_srop_single_arg(ntests, nreps, mpfr_sqrt, sqrt_flop, sr_sqrt);
  printf("all tests passed.\n");

  // Validate sr_sqrt_sr_only().
  printf("[INFO] Now testing sr_sqrt_rn_only()... ");
  randtest_srop_single_arg(ntests, nreps, mpfr_sqrt, sqrt_flop,
                           sr_sqrt_rn_only);
  printf("all tests passed.\n");


  /* Test some corner cases. */

  // Validate overflow cases of sr_add_rn_only
  printf("[INFO] Now testing sr_add_rn_only special cases... ");
  special_tests_sradd_rn_only(nreps, sr_add_rn_only);
  printf("all tests passed.\n");

  // Validate special cases of sr_add
  printf("[INFO] Now testing sr_add special cases... ");
  special_tests_sradd(nreps, sr_add);
  printf("all tests passed.\n");

  // Validate special cases of sr_mul_fma_rn_only
  printf("[INFO] Now testing sr_mul_rn_only special cases... ");
  special_cases_srmul_rn_only(nreps, sr_mul_fma_rn_only);
  printf("all tests passed.\n");

  // Validate special cases of sr_mul_fma
  printf("[INFO] Now testing sr_mul special cases... ");
  special_cases_srmul(nreps, sr_mul_fma);
  printf("all tests passed.\n");

  // Terminate.
  free(counters);
  mpfr_clears(xh, xdmp, work, (mpfr_ptr) 0);
  mpfr_free_cache ();

  return 0;
}
