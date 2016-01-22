#ifndef DWP_STATS_H
#define DWP_STATS_H

#include "filzbach.h"

//double negative_binomial_density(int, double, double, int);

double gamma(double x);
double gamma_function(double n);
double beta(double x, double y);
double factorial(int);
double gamma_density(double, double, double);
double spearman_rank(double[], double[], int);
void testing_neg_binom(void);
void test_spearman(void);
void testing_logit(void);

void initialize_stat(); // initialize draw algorithms so that they produce deterministic sequences

double gamma_draw(double, double);
double gamma_draw_shape(int, double);
double gamma_draw_mean(int, double);
double gamma_draw_stdev(int, double);
double beta_draw_ab(int, int);
#if !defined(drand)
double drand();
#endif

#endif // DWP_STATS_H
