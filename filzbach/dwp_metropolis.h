#ifndef DWP_METROPOLIS_H
#define DWP_METROPOLIS_H

#include "filzbach.h"

extern int get_currentchain();

/* iters before delayed parameters begin to be altered */
const int DELAY = 50000; 

/*
void metr_setup_param(int i, char name[], double lb, double ub, double val, int type, int fixed, int dly, int dsply);
double parameter_getvalue(int mm);
double parameter_getvalue(char name[]);
int parameter_ifchanged(char name[]);
int parameter_ifchanged(int ii);
int parameter_ifchanged(char name[], int ii);
int parameter_returnindex(char name[]);
int parameter_returnindex(char name[], int number);
void parameter_addprior(int mm, double mean, double sdev);
void parameter_copy_to_chains(void);
//TODO: remove void lnlike(void);
void lnlike_priors(void);
void calc_posteriors(void);
void calc_MLE_intervals(void);
void final_metro_output(FILE* mmfile, FILE* mmfile2);
void params_draw_random_vector(void);
void calc_IC();
void calc_convergence();

double get_metr_ltotnew(void);
long int get_metr_number_ok(void);
*/
#endif // DWP_METROPOLIS_H