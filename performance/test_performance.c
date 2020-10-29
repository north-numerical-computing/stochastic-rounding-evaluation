#include <time.h>
#include <string.h>

#include "stochrndhp.h"
#include "stochrnddw.h"

/* Structure that holds timing data for a single algorithm. */
struct timingdata {
  char method [50];
  double min;
  double max;
  double mean;
  double median;
  double deviation;
};

/* Compare two doubles, function pointer used by qsort. */
int cmpfun(const void *x, const void *y) {
  if (*(double *)x < *(double *)y)
    return -1;
  else if (*(double *)x > *(double *)y)
    return 1;
  else
    return 0;
}

/* Time new methods or stochastic rounding. */
void time_sr(char *methodname, const size_t ntests, const size_t nreps,
             struct timingdata *data,
             double (* sr_fun)(double, double)) {
  size_t i, j;
  double a, b;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX + ldexp(1., -1022);
    b = (double)rand()/RAND_MAX + ldexp(1., -1022);
    start = clock();
    for (j=0; j<nreps; j++)
      sr_fun(a, b);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i] < data->min ? timings[i] : data->min;
    data->max = timings[i] > data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  qsort(timings, ntests, sizeof(*timings), cmpfun);
  data->median = timings[ntests/2];
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
  free(timings);
}

/* Time methods for stochastic rounding based on MPFR. */
void time_sr_mpfr(char *methodname, const size_t ntests, const size_t nreps,
                  struct timingdata *data,
                  int (* mpfr_fun)(mpfr_ptr,
                                   mpfr_srcptr,
                                   mpfr_srcptr,
                                   mpfr_rnd_t),
                  mpfr_t work1, mpfr_t work2, mpfr_t work3) {
  size_t i, j;
  double a, b;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = 10202 * ((double)rand()/RAND_MAX + ldexp(1., -1022));
    b = (double)rand()/RAND_MAX + ldexp(1., -1022);
    start = clock();
    for (j=0; j<nreps; j++)
      sr_op_mpfr_noalloc(a, b, mpfr_fun, work1, work2, work3);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i]<data->min ? timings[i] : data->min;
    data->max = timings[i]>data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  qsort(timings, ntests, sizeof(*timings), cmpfun);
  data->median = timings[ntests/2];
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
  free(timings);
}

/* Time MPFR functions for reference. */
void time_mpfr(char  * methodname, const size_t ntests, const size_t nreps,
               struct timingdata *data,
               int (* mpfr_fun)(mpfr_ptr,
                                mpfr_srcptr,
                                mpfr_srcptr,
                                mpfr_rnd_t),
               mpfr_t work1, mpfr_t work2) {
  size_t i, j;
  double a, b;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX + ldexp(1., -1022);
    b = (double)rand()/RAND_MAX + ldexp(1., -1022);
    mpfr_set_d(work1, a, MPFR_RNDN);
    mpfr_set_d(work2, b, MPFR_RNDN);
    start = clock();
    for (j=0; j<nreps; j++)
      mpfr_fun(work1, work1, work2, MPFR_RNDN);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i]<data->min ? timings[i] : data->min;
    data->max = timings[i]>data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  qsort(timings, ntests, sizeof(*timings), cmpfun);
  data->median = timings[ntests/2];
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
  free(timings);
}

/* Time new methods for stochastic rounding. */
void time_sr_single_arg(char *methodname, const size_t ntests,
                        const size_t nreps,
                        struct timingdata *data,
                        double (* sr_fun)(double)) {
  size_t i, j;
  double a;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX + ldexp(1., -1022);
    start = clock();
    for (j=0; j<nreps; j++)
      sr_fun(a);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i] < data->min ? timings[i] : data->min;
    data->max = timings[i] > data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  qsort(timings, ntests, sizeof(*timings), cmpfun);
  data->median = timings[ntests/2];
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
}

/* Time methods for stochastic rounding based on MPFR. */
void time_sr_mpfr_single_arg(char *methodname,
                             const size_t ntests, const size_t nreps,
                             struct timingdata *data,
                             int (* mpfr_fun)(mpfr_ptr,
                                              mpfr_srcptr,
                                              mpfr_rnd_t),
                             mpfr_t work1, mpfr_t work2, mpfr_t work3) {
  size_t i, j;
  double a;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX + ldexp(1., -1022);
    start = clock();
    for (j=0; j<nreps; j++)
      sr_op_mpfr_noalloc_single_arg(a, mpfr_fun, work1, work2, work3);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i]<data->min ? timings[i] : data->min;
    data->max = timings[i]>data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  qsort(timings, ntests, sizeof(*timings), cmpfun);
  data->median = timings[ntests/2];
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
  free(timings);
}

/* Time MPFR functions for reference. */
void time_mpfr_single_arg(char *methodname,
                          const size_t ntests, const size_t nreps,
                          struct timingdata *data,
                          int (* mpfr_fun)(mpfr_ptr,
                                           mpfr_srcptr,
                                           mpfr_rnd_t),
                          mpfr_t work1) {
  size_t i, j;
  double a;
  clock_t start;
  strcpy(data->method, methodname);
  data->min = INFINITY;
  data->max = -INFINITY;
  data->mean = 0;
  double *timings = malloc(ntests * sizeof(double));
  for (i=0; i<ntests; i++) {
    a = (double)rand()/RAND_MAX + ldexp(1., -1022);
    mpfr_set_d(work1, a, MPFR_RNDN);
    start = clock();
    for (j=0; j<nreps; j++)
      mpfr_fun(work1, work1, MPFR_RNDN);
    timings[i] = (double)(CLOCKS_PER_SEC * nreps) / (clock() - start);
    data->min = timings[i]<data->min ? timings[i] : data->min;
    data->max = timings[i]>data->max ? timings[i] : data->max;
    data->mean += timings[i];
  }
  data->mean = data->mean / ntests;
  data->deviation = 0.;
  for (i=0; i<ntests; i++)
    data->deviation += pow(timings[i] - data->mean, 2);
  data->deviation = sqrt(data->deviation/(ntests-1));
  free(timings);
}

void print_annextable_row(FILE *file, struct timingdata *data) {
  fprintf(file, "| %13s | %8.2f | %8.2f | %8.2f | %8.2f |\n",
          data->method,
          data->min/1e6, data->max/1e6,
          data->mean/1e6, data->deviation/1e6);
}
void print_annextable_rule(FILE *file) {
  fprintf(file,
          "+---------------+----------+----------+----------+----------+\n");
}

int main () {

  // Initialize PRNGs.
  const unsigned long seed = 1;
  srand(seed);

  // Initialize timing loops parameters.
  size_t ntests = 100;
  size_t nreps = 10000000; // 1e7
  size_t nalgs = 20;

  struct timingdata data [nalgs];

  printf("Running %ld tests, %ld repetitions each.\n", ntests, nreps);

  // New algorithms.
  printf("[INFO] Running new algorithms... ");

  time_sr("sr_add", ntests, nreps, &data[3], sr_add);
  time_sr("sr_add_rn_only", ntests, nreps, &data[4], sr_add_rn_only);

  time_sr("sr_mul_fma", ntests, nreps, &data[8], sr_mul_fma);
  time_sr("sr_mul_fma_rn_only", ntests, nreps, &data[9], sr_mul_fma_rn_only);

  time_sr("sr_div", ntests, nreps, &data[13], sr_div);
  time_sr("sr_div_rn_only", ntests, nreps, &data[14], sr_div_rn_only);

  time_sr_single_arg("sr_sqrt", ntests, nreps, &data[18], sr_sqrt);
  time_sr_single_arg("sr_sqrt_rn_only", ntests, nreps, &data[19],
                     sr_sqrt_rn_only);

  printf("done.\n");

  // MPFR-based algorithms.
  size_t i;
  size_t p [3] = {61, 88, 113};
  mpfr_t work1, work2, work3;
  char* algname = malloc(50 * sizeof(char));
  for (i=0; i<3; i++) {
    printf("[INFO] Running algorithms in %lu bits... ", p[i]);
    mpfr_set_default_prec(p[i]);
    mpfr_inits2(p[i], work1, work2, work3, (mpfr_ptr) 0);
    sprintf(algname, "sr_add_mpfr(%lu)", p[i]);
    time_sr_mpfr(algname, ntests, nreps, &data[i],
                 mpfr_add, work1, work2, work3);
    sprintf(algname, "sr_mul_mpfr(%lu)", p[i]);
    time_sr_mpfr(algname, ntests, nreps, &data[5+i],
                 mpfr_mul, work1, work2, work3);
    sprintf(algname, "sr_div_mpfr(%lu)", p[i]);
    time_sr_mpfr(algname, ntests, nreps, &data[10+i],
                 mpfr_div, work1, work2, work3);
    sprintf(algname, "sr_sqrt_mpfr(%lu)", p[i]);
    time_sr_mpfr_single_arg(algname, ntests, nreps, &data[15+i],
                            mpfr_sqrt, work1, work2, work3);
    mpfr_clears(work1, work2, work3, (mpfr_ptr) 0);
    printf("done.\n");
  }
  free(algname);

  // Print to stdout LaTeX-ready table.
  char *terminator = malloc(10 * sizeof(*terminator));
  FILE * f = stdout;

  size_t reference;
  fprintf(f, "min       & ");
  for (i=0; i<nalgs; i++) {
    terminator = i==5||i==10||i==15? "&& " : i<nalgs-1? "& " : "";
    fprintf(f, "%s", terminator);
    fprintf(f, "%4.2f\\hphantom{$\\times$}", data[i].min/1e6);
  }
  fprintf(f, "\\\\\n");

  fprintf(f, "max       & ");
  for (i=0; i<nalgs; i++) {
    terminator = i==5||i==10||i==15? "&& " : i<nalgs-1? "& " : "";
    fprintf(f, "%s", terminator);
    fprintf(f, "%4.2f\\hphantom{$\\times$}", data[i].max/1e6);
  }
  fprintf(f, "\\\\\n");

  fprintf(f, "mean      & ");
  for (i=0; i<nalgs; i++) {
    terminator = i==5||i==10||i==15? "&& " : i<nalgs-1? "& " : "";
    fprintf(f, "%s", terminator);
    fprintf(f, "%4.2f\\hphantom{$\\times$}", data[i].mean/1e6);
  }
  fprintf(f, "\\\\\n");

  fprintf(f, "\\rowstyle{}\\hspace{2pt}$\\hookrightarrow$\\hspace{2pt}\\textit{speedup} & ");
  for (i=0; i<nalgs; i++) {
    reference = i<5 ? 2 : (i<10 ? 7 : (i<15 ? 12 : 17));
    terminator = i==5||i==10||i==15? "&& " : i<nalgs-1? "& " : "";
    fprintf(f, "%s", terminator);
    fprintf(f, "%4.2f$\\times$", data[i].mean/data[reference].mean);
  }
  fprintf(f, "\\\\\n");

  fprintf(f, "deviation & ");
  for (i=0; i<nalgs; i++) {
    terminator = i==5||i==10||i==15? "&& " : i<nalgs-1? "& " : "";
    fprintf(f, "%s", terminator);
    fprintf(f, "%4.2f\\hphantom{$\\times$}", data[i].deviation/1e6);
  }
  fprintf(f, "\\\\\n");

  // Print to stdout pretty table.
  // AR = use all rounding modes
  // RN = use only round-to-nearest
  fprintf(f, "\n");
  fprintf(f, "+-----------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+\n");
  fprintf(f, "| Algorithm |");
  fprintf(f, "   +61 |   +88 |  +113 |   +AR |   +RN |");
  fprintf(f, "   ×61 |   ×88 |  ×113 |   ×AR |   ×RN |");
  fprintf(f, "   ÷61 |   ÷88 |  ÷113 |   ÷AR |   ÷RN |");
  fprintf(f, "   √61 |   √88 |  √113 |   √AR |   √RN |\n");
  fprintf(f, "+-----------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+\n");
  fprintf(f, "| min       | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    fprintf(f, "%5.1f", data[i].min/1e6);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| max       | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    fprintf(f, "%5.1f", data[i].max/1e6);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| mean      | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    fprintf(f, "%5.1f", data[i].mean/1e6);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| ↳speedup  | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    reference = i<5? 2: (i<10? 7: 12);
    fprintf(f, "%5.1f", data[i].mean/data[reference].mean);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| deviation | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    fprintf(f, "%5.1f", data[i].deviation/1e6);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| median    | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    fprintf(f, "%5.1f", data[i].median/1e6);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "| ↳speedup  | ");
  for (i=0; i<nalgs; i++) {
    terminator = i<nalgs-1? " | " : " |\n";
    reference = i<5? 2: (i<10? 7: 12);
    fprintf(f, "%5.1f", data[i].median/data[reference].median);
    fprintf(f, "%s", terminator);
  }
  fprintf(f, "+-----------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+");
  fprintf(f, "-------+-------+-------+-------+-------+\n\n");

  // Test performance of MPFR operations, for comparison with the data in:
  // Fabiano, Muller, & Picot, "Algorithms for Triple-Word Arithmetic,"
  // IEEE Trans. Comp., 68(11):1573-1583, 2019. DOI: 10.1109/TC.2019.2918451
  size_t pp = 159; // Precision bits in quad precision.
  struct timingdata adddata = {};
  mpfr_set_default_prec(pp);
  print_annextable_rule(stdout);
  printf("|        method |      min |      max |     mean | std. dev.|\n");
  print_annextable_rule(stdout);
  mpfr_inits2(pp, work1, work2, (mpfr_ptr) 0);
  time_mpfr("mpfr_add", ntests, nreps, &adddata,
            mpfr_add, work1, work2);
  print_annextable_row(stdout, &adddata);
  time_mpfr("mpfr_sub", ntests, nreps, &adddata,
            mpfr_sub, work1, work2);
  print_annextable_row(stdout, &adddata);
  time_mpfr("mpfr_mul", ntests, nreps, &adddata,
            mpfr_mul, work1, work2);
  print_annextable_row(stdout,  &adddata);
  time_mpfr("mpfr_div", ntests, nreps, &adddata,
            mpfr_div, work1, work2);
  print_annextable_row(stdout, &adddata);
  time_mpfr_single_arg("mpfr_sqrt", ntests, nreps, &adddata,
                       mpfr_sqrt, work1);
  print_annextable_row(stdout, &adddata);
  print_annextable_rule(stdout);
  mpfr_clears(work1, work2, (mpfr_ptr) 0);

  // Terminate.
  mpfr_free_cache ();

  return 0;
}
