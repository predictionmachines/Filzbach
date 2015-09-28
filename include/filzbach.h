#ifndef FILZBACH_H
#define FILZBACH_H

#ifndef FILZBACH_API
#ifdef FILZBACH_DLL
 #define FILZBACH_API __declspec(dllexport)
 #define EXTERN_FILZBACH_API extern __declspec(dllexport)
#elif defined(FILZBACH_LIB)
 #define FILZBACH_API extern
 #define EXTERN_FILZBACH_API extern
#else
 #define FILZBACH_API __declspec(dllimport)
 #define EXTERN_FILZBACH_API extern __declspec(dllimport)
#endif
#endif

#include <cstddef>

FILZBACH_API void table_writevalue(const char table_name[], const char column_name[], int row_number, double value_to_write);
FILZBACH_API void table_writevalue(const char table_name[], int column_number, int row_number, double value_to_write);
FILZBACH_API void table_writevalue_multichain(const char table_name[], const char column_name[], int row_number, double value_to_write);
FILZBACH_API void table_writevalue_multichain(const char table_name[], int column_number, int row_number, double value_to_write);
FILZBACH_API int table_create(const char table_name[]);
FILZBACH_API int table_addcolumn(const char table_name[], const char column_name[]);
FILZBACH_API int table_addrows(const char table_name[], int number_of_rows_to_add);
FILZBACH_API int table_removerows(const char table_name[], int number_of_rows_to_delete);
FILZBACH_API int table_getcolumnID(int  table_numeric_id, const char column_name[]);
FILZBACH_API void table_setnumrows(const char table_name[],int number_of_rows);
FILZBACH_API void table_read(const char table_name[] , const char path_and_file_name[], int number_of_columns);
FILZBACH_API int table_numrows(const char table_name[]);
FILZBACH_API int table_numrows(int table_numeric_id);
FILZBACH_API int table_numcolumns(const char table_name[]);
FILZBACH_API int table_numcolumns(int table_numeric_id);
FILZBACH_API double table_getvalue(const char table_name[], const char column_name[], int row_number);
FILZBACH_API double table_getvalue(const char table_name[], int column_number, int row_number);
FILZBACH_API void table_output(const char table_name[], const char path_and_file_name[]);
FILZBACH_API int table_getID(const char table_name[]);
FILZBACH_API void table_print(const char table_name[]);
FILZBACH_API int table_delete(const char table_name[]);
FILZBACH_API int table_delete(int table_numeric_id);
FILZBACH_API int table_deletecolumn(const char table_name[], const char column_name[]);
FILZBACH_API int table_deletecolumn(const char table_name[], int column_number);
FILZBACH_API int table_deletecolumn(int table_numeric_id, const char column_name[]);
FILZBACH_API int table_deletecolumn(int table_numeric_id, int column_number);
FILZBACH_API double table_getcolumnmax(const char table_name[], const char column_name[]);
FILZBACH_API double table_getcolumnmax(const char table_name[], int column_number);
FILZBACH_API double table_getcolumnmin(const char table_name[], const char column_name[]);
FILZBACH_API double table_getcolumnmin(const char table_name[], int column_number);
FILZBACH_API int table_getvaluecount(const char table_name[], const char column_name[]);
FILZBACH_API int table_getvaluecount(const char table_name[], int column_number);

EXTERN_FILZBACH_API void (*pfn_likelihood) (void);
EXTERN_FILZBACH_API double (*pfn_parallel_likelihood)(int, int, double *ltot, long int *numok);

FILZBACH_API void set_parallel_likelihood(bool, int, int);

FILZBACH_API void set_output_options(int console, int summary, int params, int bayes, int mle, int chain);
FILZBACH_API void initialize_filzbach(void);
FILZBACH_API void set_chains(int);
FILZBACH_API void set_thinning(int);
FILZBACH_API void runmcmc(int tburnin, int teststeps, int tburnin2, int tmleexp);
FILZBACH_API void parameter_create_vector(const char name[], double lb, double ub, double val, int type, int fixed, int dsply, int number);
FILZBACH_API int parameter_create(const char name[], double lb, double ub, double val, int type, int fixed, int dsply);
FILZBACH_API void parameter_showall(void);
FILZBACH_API double parameter_getvalue(const char name[]);
FILZBACH_API double parameter_getvalue(const char name[], int number);
FILZBACH_API void parameter_setvalue(const char name[], double val);
FILZBACH_API void parameter_setvalue(int mm, double val);
FILZBACH_API void parameter_setvalue(const char name[], int index, double val);
FILZBACH_API double cv(const char[]);
FILZBACH_API double cv(int);
FILZBACH_API double cv(const char[], int);
FILZBACH_API double l95(const char name[]);
FILZBACH_API double l95(const char name[], int number);
FILZBACH_API double l68(const char name[]);
FILZBACH_API double l68(const char name[], int number);
FILZBACH_API double postmedian(const char name[], int number);
FILZBACH_API double postmedian(const char name[]);
FILZBACH_API double postmean(const char name[]);
FILZBACH_API double postmean(const char name[], int number);
FILZBACH_API double u68(const char name[]);
FILZBACH_API double u68(const char name[], int number);
FILZBACH_API double u95(const char name[]);
FILZBACH_API double u95(const char name[], int number);
FILZBACH_API double pl_l95(const char name[]);
FILZBACH_API double pl_l95(const char name[], int number);
FILZBACH_API double pl_l68(const char name[]);
FILZBACH_API double pl_l68(const char name[], int number);
FILZBACH_API double mle(const char name[], int number);
FILZBACH_API double mle(const char name[]);
FILZBACH_API double pl_u68(const char name[]);
FILZBACH_API double pl_u68(const char name[], int number);
FILZBACH_API double pl_u95(const char name[]);
FILZBACH_API double pl_u95(const char name[], int number);
FILZBACH_API int parameter_ifchanged(const char name[]);
FILZBACH_API int parameter_ifchanged(int ii);
FILZBACH_API int parameter_ifchanged(const char name[], int ii);
FILZBACH_API int parameter_returnindex(const char name[]);
FILZBACH_API int parameter_returnindex(const char name[], int number);
FILZBACH_API void parameter_addprior(const char name[], double mean, double sdev);
FILZBACH_API void parameter_addprior_vector(const char name[], double mean, double sdev);
FILZBACH_API void parameter_normalize_vector(const char name[], double mean, double sdev);
FILZBACH_API void name_analysis(const char name[]);
FILZBACH_API void params_set_to_posterior_mean(void);
FILZBACH_API void set_metr_ltotnew(double);
FILZBACH_API void set_metr_number_ok(long int);
FILZBACH_API void inc_metr_ltotnew(double);
FILZBACH_API void inc_metr_number_ok(long int);
FILZBACH_API void get_filzbach_path(char dest[], size_t bufSize);
FILZBACH_API void force_accept(void);
FILZBACH_API void force_reject(void);
FILZBACH_API void params_set_to_old(void);
FILZBACH_API void params_draw_random_vector(void);
FILZBACH_API int params_from_bayes_list(int index);
FILZBACH_API int params_from_bayes_list(int chain, int index);
FILZBACH_API int inburnin();
FILZBACH_API void parameter_delay(const char name[], int iters);
FILZBACH_API void parameter_fix(const char name[]);
FILZBACH_API void parameter_fix(const char name[], int number);
FILZBACH_API void parameter_fix_value(const char name[], double val);
FILZBACH_API void parameter_fix_value(const char name[], int number, double val);
FILZBACH_API void parameter_unfix(const char name[]);
FILZBACH_API void parameter_unfix(const char name[], int number);


FILZBACH_API double normal_density(double observation, double mean, double sdev);
FILZBACH_API double normal_draw(double mean, double sdev);

FILZBACH_API double gaussian(double x, double optimum, double width_parameter);

FILZBACH_API double poisson_density(int count, double mean);
FILZBACH_API int poisson_draw(double mean);

FILZBACH_API double gamma_density(double observation, double mean, double sdev);
// No corresponding draw function available 

FILZBACH_API double gamma_density_scale(double observation, double mean, double scale_parameter);
// No corresponding draw function available 

 FILZBACH_API double gamma_density_shape(double observation, double mean, double shape_parameter);
// No corresponding draw function available 

FILZBACH_API double binomial_density(int heads, int number_trials, double prob_head_for_each_trial);
FILZBACH_API double binomial_draw(int number_trials, double prob_head_for_each_trial);

FILZBACH_API double negative_binomial_density_shape(int count, double mean, double shape_parameter);
// No corresponding draw function available 

FILZBACH_API double negative_binomial_density(int count, double mean, double sdev);
FILZBACH_API double negative_binomial_draw(double mean, double sdev);

FILZBACH_API double negative_binomial_density_p(int k, int r, double p);
FILZBACH_API double negative_binomial_draw_p(double, double);

FILZBACH_API double beta_density(double observation, double mean, double sdev);
// No corresponding draw function available -- this function is bad and will be removed soon

FILZBACH_API double beta_density_mean_fvar(double observation, double mean, double fvar);
// No corresponding draw function available -- this is the new version that should work

FILZBACH_API double beta_density_ab(double observation, double parameter_a, double parameter_b);
// No corresponding draw function available -- bad function, to be removed

FILZBACH_API double beta_density_alpha_beta(double observation, double alpha, double beta);
// No corresponding draw function available -- good version

FILZBACH_API double beta_draw_mean_fvar(double mean, double fvar);
// No corresponding draw function available -- good version

FILZBACH_API double beta_draw_alpha_beta(double alpha, double beta);
// No corresponding draw function available -- good version

FILZBACH_API double beta_draw_mean_alpha(double mean, double alpha);
// No corresponding draw function available -- good version

FILZBACH_API double beta_density_mean_alpha(double observation, double mean, double alpha);
// No corresponding draw function available 

FILZBACH_API double exponential_density(double observation, double mean);
FILZBACH_API double exponential_draw(double mean);

FILZBACH_API double exponential_rate(double observation, double rate_parameter);
FILZBACH_API double exponential_rate_draw(double rate_parameter);

//FILZBACH_API double lognormal_density(double observation, double mean, double sdev);
//FILZBACH_API double lognormal_draw(double mean, double sdev);

FILZBACH_API double lognormal_density(double x, double mode, double stdev_in_logspace);
FILZBACH_API double lognormal_draw(double mode, double stdev_in_logspace);

FILZBACH_API double logistic(double k);
FILZBACH_API double random(double lower_bound, double upper_bound);
FILZBACH_API int random_integer(int lower_bound, int upper_bound);
FILZBACH_API double large(double value[], int lid, int uid, int num);
FILZBACH_API double large(float value[], int lid, int uid, int num);
FILZBACH_API double large(int value[], int lid, int uid, int num);

#endif // FILZBACH_H
