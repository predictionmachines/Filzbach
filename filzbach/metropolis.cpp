#include "stdafx.h"
#include <cassert>
#include "dwp_stats.h"
#include "dwp_metropolis.h"
#include "params.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#define _chdir chdir
inline int _mkdir(const char* dirname) { mkdir(dirname, 0777); }
#define strcpy_s(dest,bufsize,source) strncpy(dest,source,bufsize)
#define strcat_s(dest,bufsize,source) strncat(dest,source,bufsize)
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <ctime>

#include <vector>
#include <algorithm>
#include <functional>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Managed Likelihood Callback function.
// This function is used to set the callback to ptrlnlike
// from managed code.
//
// The corresponding import in C# would be:
//
// [DllImport("Filzbach.dll", CallingConvention = CallingConvention.Cdecl, EntryPoint = "?likelihood_callback@@YAXP6AXXZ@Z")]
// public static extern void SetLikelihood(LikelihoodCallback callback);
//
// Then create a delegate and the corresponding callback method.
//
// public delegate void LikelihoodCallback();
// public void Likelihood()
// {
//   ... your codes goes here
// }
//
// To set the callback frm managed code use the following lines:
// The fist and last line is to pinn the delegate
// so the GC does not clean up the delegate
//
// LikelihoodCallback callback = this.Likelihood;
// SetLikelihood(callback);
// ... your code goes here
// callback = null;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void(*pfn_likelihood) (void);
double(*pfn_parallel_likelihood)(int, int, double *ltot, long int *numok);

typedef void(*PFN_LIKELIHOOD_CALLBACK)(void);
FILZBACH_API void likelihood_callback(PFN_LIKELIHOOD_CALLBACK callback);
extern void likelihood_callback(PFN_LIKELIHOOD_CALLBACK callback)
{
	pfn_likelihood = (void(*)(void))callback;
}

typedef double(*PFN_PARALLEL_LIKELIHOOD_CALLBACK)(int start, int end, double *ltot, long int *numok);
EXTERN_FILZBACH_API void parallel_likelihood_callback(PFN_PARALLEL_LIKELIHOOD_CALLBACK callback)
{
	pfn_parallel_likelihood = (double(*)(int, int, double *, long int *))callback;
}


/* number of steps between screen output */
int console_out_step = 1000;
// steps for computing bayesian values used within calc_posteriors
int bayes_step = 100;

//execute likelihood in parallel
bool parallel_mode = false;
// iterations for parallel likelihood
int likelihood_steps = 0;
int likelihood_lbound = 0;
int likelihood_ubound = 0;

// number of chains
int chaincount = 1;
// tables that hold the bayesian value
double ***bayestable = NULL;
// length of bayestable (that but containing only the bayesian samples part)
int bayestable_length = 0;
// length of space to hold bayestable info (this also contains samples from extra phases to do lhood profile method from phases 3, 4)
int bayestable_length_extras = 0;
// position when processing bayesian values for chains
int* bayespos = NULL;
int* bayes_chain_length = NULL;
// indicates the currently processed chain, -1 if no chain is processed, thread private variable
int currentchain = -1;
#pragma omp threadprivate(currentchain)
int get_currentchain() { return currentchain; }

int metr_k = 1;
#pragma omp threadprivate(metr_k)
int bestchainID = 0; // this value is overwritten after the MCMC sampling has finished
// 
int *numalt = NULL;

double *aic = NULL;
double *bic = NULL;
double *dic = NULL;
double base_aic;
double base_bic;
double baseDic;
double convergence_RhatMean; // each param gets an Rhat value -- this is the mean of these

int base_freeparams;
int *freeparams = NULL;

double base_ltotmax;
double *ltotmax = NULL;

double base_ltotnew;
double *ltotnew = NULL;

double base_ltotold;
double *ltotold = NULL;

double base_ptotnew;
double *ptotnew = NULL;

double base_ptotold;
double *ptotold = NULL;

// chain specific 
int *accept_forced = NULL;
int *reject_forced = NULL;

// counter variables
double *pbacc = NULL;
double *pbcount = NULL;

/**********************************************
 * Output options
 **********************************************/

 // analysis name
char metr_tag[256] = "Default";
bool out_summary = true;
bool out_params = true;
bool out_bayes = true;
bool out_mle = true;
bool out_chain = false;

#define out_summary_name "_MCMC_final_out.txt"
#define out_params_name "_MCMC_params_out.txt"
#define out_bayes_name "_MCMC_bayes_list.txt"
#define out_mle_name "_MCMC_MLE_list.txt"
#define out_chain_name "_MCMC_metr_out.txt"

#define output_option(name) if (name >= 0) out_##name = name>0;
void set_output_options(int console, int summary, int params, int bayes, int mle, int chain)
{
	if (console >= 0) console_out_step = console;
	output_option(summary)
		output_option(params)
		output_option(bayes)
		output_option(mle)
		output_option(chain)
}

void get_filzbach_path(char dest[], size_t bufSize)
{
	strcpy_s(dest, bufSize, ws_dir);
	strcat_s(dest, bufSize, "/");
	strcat_s(dest, bufSize, metr_tag);
}

FILE* workspace_fopen(const char* name, const char* options)
{
	char outbuf[512];
	get_filzbach_path(outbuf, 512);
	strcat(outbuf, name);
#pragma warning(push)
#pragma warning(disable:4996)
	FILE* result = fopen(outbuf, options);
#pragma warning(pop)
	CHECK(result,"couldn't open output file")
	return result;
}

/**********************************************
 * metropolis algorithm specific variables
 **********************************************/

int burnin;
int mleexp;
int burnin2;
int eststeps;

long int base_metr_number_ok;
long int *metr_number_ok = NULL;


/*
 * Sets the parallel execution mode for the parallel likelihood.
 */
void set_parallel_likelihood(bool parallel, int lbound, int ubound)
{
	parallel_mode = parallel;
	likelihood_lbound = lbound;
	likelihood_ubound = ubound;
}

/*
*
*/
double get_metr_ltotnew()
{
	return currentchain < 0 ? base_ltotnew : ltotnew[currentchain];
}

/*
*
*/
void set_metr_ltotnew(double value)
{
	if (currentchain < 0)
	{
		base_ltotnew = value;
	}
	else
	{
		ltotnew[currentchain] = value;
	}
}

/*
*
*/
void inc_metr_ltotnew(double value)
{
	if (currentchain < 0)
	{
		base_ltotnew += value;
	}
	else
	{
		ltotnew[currentchain] += value;
	}
}

/*
*
*/
long int get_metr_number_ok()
{
	return currentchain < 0 ? base_metr_number_ok : metr_number_ok[currentchain];
}

/*
*
*/
void set_metr_number_ok(long int value)
{
	if (currentchain < 0)
		base_metr_number_ok = value;
	else
		metr_number_ok[currentchain] = value;
}

void inc_metr_number_ok(long int value)
{
	if (currentchain < 0)
		base_metr_number_ok += value;
	else
		metr_number_ok[currentchain] += value;
}




/*
* For keeping track of bayesian posterior averages.
* Initializes the data structures used for the bayesian averages.
*
* teststeps - The number of teststeps for the metropolis algorithm.
*             Used to calculate the number of rows required to store
+             the bayes values.
*/
void init_bayestable(int space1, int space2)
{
	// deallocate previously allocated space
	static int bayestable_chaincount = 0;
	if ((bayestable_chaincount > 0) && (bayestable != NULL))
	{
		for (int c = 0; c < bayestable_chaincount; c++)
		{
			for (int i = 0; i < bayestable_length_extras; i++) delete[] bayestable[c][i];
			delete[] bayestable[c];
		}
		delete[] bayestable;
	}
	// allocate new tables
	bayestable_length = space1;
	bayestable_length_extras = space1 + space2;

	// init bayes tables for chains
	bayestable = new double **[chaincount];
	CHECK(bayestable != NULL, "Memory allocation failed.");
	bayestable_chaincount = chaincount;
	// create one table per chain
	for (int c = 0; c < chaincount; c++)
	{
		bayestable[c] = new double *[bayestable_length_extras]; // therefore the table is now big enough to also include the extras
		CHECK(bayestable[c] != NULL, "Memory allocation failed.");

		// initialize each entry
		for (int i = 0; i < bayestable_length_extras; i++)
		{
			// +3 for the first column + 1 for metr_k, ltotnew + ptotnew, and ltotnew			
			bayestable[c][i] = new double[paramcount + 3];
			CHECK(bayestable[c][i] != NULL, "Memory allocation failed.");

		}
	}


	// init base positions for chains
	delete[] bayespos;
	bayespos = new int[chaincount];
	CHECK(bayespos != NULL, "Memory allocation failed.");
	for (int c = 0; c < chaincount; c++)
	{
		bayespos[c] = 0;
	}

	// allocate lengths
	delete[] bayes_chain_length;
	bayes_chain_length = new int[chaincount];
	CHECK(bayes_chain_length != NULL, "Memory allocation failed.");
	for (int c = 0; c < chaincount; c++)
	{
		bayes_chain_length[c] = 0;
	}

}

/*
* For keeping track of bayesian posterior averages.
* Fills the datastructure for the bayesian posterior averages.
* Each chain requires its own table.
*/
void fill_bayestable(int chain, long int iter)
{
	ASSERT(chain >= 0 && chain < chaincount, "Index out of range")
		int pos;

	// chain-specific bayestable
	pos = bayespos[chain];

	if (pos < bayestable_length_extras)
	{

		bayestable[chain][pos][0] = iter;
		bayestable[chain][pos][1] = ltotnew[chain] + ptotnew[chain];
		bayestable[chain][pos][2] = ltotnew[chain];
		param* mpara = current_para(chain);
		for (int i = 3; i < paramcount + 3; i++)
		{
			bayestable[chain][pos][i] = mpara[i - 3].value;
		}

		bayespos[chain]++; // so after last time bayespos will equal e.g., 51, if the samples are for positions 0 - 50
	}
	else
	{
		printf("\n fill chain-specific bayestable exceeding table length, iteration %ld, iteration %d pos is %d \n", iter, chain, bayespos[chain]);
		// system("pause");
	}
}


/*
* Set the number of chains.
*/
void set_chains(int count)
{
	if (count < 1)
		return;

	chaincount = count;
}

void set_thinning(int number)
{
	if (number < 1)
		return;

	bayes_step = number;
}

void initialize_filzbach()
{
	initialize_stat();
	paramcount = 0;
	chaincount = 1;

	if (_chdir(ws_dir))
		_mkdir(ws_dir);
	else
		_chdir("..");
}

/*
 * Initializes all l and p varialbes.
 */
void init_tots()
{
	base_ltotmax = 0.0;
	delete[] ltotmax;
	ltotmax = new double[chaincount];
	CHECK(ltotmax != NULL, "Memory allocation failed.");

	base_ltotnew = 0.0;
	delete[] ltotnew;
	ltotnew = new double[chaincount];
	CHECK(ltotnew != NULL, "Memory allocation failed.");

	base_ltotold = 0.0;
	delete[] ltotold;
	ltotold = new double[chaincount];
	CHECK(ltotold != NULL, "Memory allocation failed.");

	base_ptotnew = 0.0;
	delete[] ptotnew;
	ptotnew = new double[chaincount];
	CHECK(ptotnew != NULL, "Memory allocation failed.");

	base_ptotold = 0.0;
	delete[] ptotold;
	ptotold = new double[chaincount];
	CHECK(ptotold != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; c++)
	{
		ltotmax[c] = 0.0;
		ltotnew[c] = 0.0;
		ltotold[c] = 0.0;
		ptotnew[c] = 0.0;
		ptotold[c] = 0.0;
	}
}

/*
 * Initializes counter variables only used within runmcmc.
  */
void init_counters()
{
	delete[] metr_number_ok;
	metr_number_ok = new long int[chaincount];
	CHECK(metr_number_ok != NULL, "Memory allocation failed.");

	delete[] pbacc;
	pbacc = new double[chaincount];
	CHECK(pbacc != NULL, "Memory allocation failed.");

	delete[] pbcount;
	pbcount = new double[chaincount];
	CHECK(pbcount != NULL, "Memory allocation failed.");

	base_metr_number_ok = 0;
	for (int c = 0; c < chaincount; ++c)
	{
		metr_number_ok[c] = 0;
		pbacc[c] = 0.0;
		pbcount[c] = 0.0;
	}

}

/***************************************************************************/
void name_analysis(const char name[])
{
	CHECK(name != NULL, "Analysis name cannot be null.");
	strcpy(metr_tag, name);
}



/*
* Clear alt and altt lists.
*/
void clear_altlists()
{
	for (int c = 0; c < chaincount; c++)
	{
		param* mpara = current_para(c);
		for (int i = 0; i < paramcount; i++)
		{
			mpara[i].alt = 0;
			mpara[i].altt = 0;
			mpara[i].runalt = 0;
			mpara[i].runacc = 0;
		}
	}
}

/*
* Counts the number of fixed paramters for each chain.
*/
int* countFreeParams()
{
	int *freeParamCount = new int[chaincount];
	CHECK(freeParamCount != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; c++)
	{
		freeParamCount[c] = 0;

		param* mpara = current_para(c);
		for (int i = 0; i < paramcount; i++)
		{
			if (mpara[i].fixed < 1)
				freeParamCount[c]++;
		}
	}

	return freeParamCount;
}

/*
* TODO: Add some meaningfull description here.
* TODO: Maybe rename the function once  known what talt means.
*/
double* calcTAlt(int* fixedParams)
{
	double *talt = new double[chaincount];
	CHECK(talt != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; c++)
	{
		talt[c] = 3.0 / (double)fixedParams[c];
	}

	return talt;
}

/*
* Initialise step sizes.
*/
void initStepSizes()
{
	for (int c = 0; c < chaincount; c++)
	{
		param* mpara = current_para(c);
		for (int i = 0; i < paramcount; i++)
		{
			if (mpara[i].type < 1)
			{
				mpara[i].delta = 0.50*(mpara[i].ub - mpara[i].lb);
			}

			if (mpara[i].type > 0)
			{
				mpara[i].delta = 0.50*(log(mpara[i].ub) - log(mpara[i].lb)); // new addition (Feb 14 2011) to make use of ub,lb information
			}
		}
	}
}

/*
* Set initial values for parameters at random.
*/
void initRandomValues()
{
	for (int c = 0; c < chaincount; c++)
	{
		param* mpara = current_para(c);
		for (int i = 0; i < paramcount; i++)
		{
			if ((mpara[i].fixed == 0) && (mpara[i].delay < 1))
			{
				mpara[i].value = random(mpara[i].lb, mpara[i].ub);
			}
		}
	}
}

/*
* Initialise parameter information.
*/
void init_paramerniformation()
{
	for (int c = 0; c < chaincount; c++)
	{
		param* mpara = current_para(c);
		for (int i = 0; i < paramcount; i++)
		{
			mpara[i].pold = mpara[i].value;
			mpara[i].MLE = mpara[i].value;
		}
	}
}

/*
* Initialize the free parameter sets.
*/
void init_freeparams()
{
	base_freeparams = 0;
	delete[] freeparams;
	freeparams = new int[chaincount];
	CHECK(freeparams != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; ++c)
	{
		freeparams[c] = 0;
	}
}

/*
 * Initializes the IC values.
 */
void init_ic()
{
	base_aic = 0;
	base_bic = 0;
	baseDic = 0;

	delete[] aic;
	aic = new double[chaincount];
	CHECK(aic != NULL, "Memory allocation failed.");

	delete[] bic;
	bic = new double[chaincount];
	CHECK(bic != NULL, "Memory allocation failed.");

	delete[] dic;
	dic = new double[chaincount];
	CHECK(dic != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; ++c)
	{
		aic[c] = 0.0;
		bic[c] = 0.0;
		dic[c] = 0.0;
	}
}



/***************************************************************************/
void lnlike_priors()
{
	ASSERT(currentchain >= 0, "Internal error");
	double prob1;

	ptotnew[currentchain] = 0.0;

	param* params = current_para();
	/* loops over all params and evaluates probability against priors */
	for (int mm = 0; mm < paramcount; mm++)
	{
		if (params[mm].prioryn > 0)
		{
			if (params[mm].type == 0)
			{
				prob1 = normal_density(params[mm].value, params[mm].priormean, params[mm].priorsdev);
			}
			else
			{
				prob1 = normal_density(log(params[mm].value), log(params[mm].priormean), params[mm].priorsdev);
			}
			ptotnew[currentchain] += log(prob1);
		}
	}
}


int bestchain()
{
	// choose 'best' chain -- simply using average likelihood from sampling phase

	int bestcc = 0;
	double max_av_lhood = 0;

	// loop over all chains
	for (int cc = 0; cc < chaincount; cc++)
	{
		// calc av likelihood from samples
		int n = bayes_chain_length[cc];
		double av_lhood = 0.0;
		for (int i = 0; i < n; i++)
		{
			av_lhood += bayestable[cc][i][1];
		}
		av_lhood /= (double)n;

		if (cc == 0 || av_lhood > max_av_lhood)
		{
			max_av_lhood = av_lhood;
			bestcc = cc;
		}
	}

	return bestcc;
}

/*
 * Consolidating chains and selecting the chains that should be referred as base chain after computation.
 */
void consolidateChains()
{
	bestchainID = bestchain();
	base_freeparams = freeparams[bestchainID];
	int chain = bestchainID;

	param* base_para = current_para(-1);
	param* params = current_para(chain);
	for (int p = 0; p < paramcount; ++p)
	{
		base_para[p] = params[p];
	}

	base_aic = aic[chain];
	base_bic = bic[chain];
	baseDic = dic[chain];

	base_ltotmax = ltotmax[chain];
	base_ltotold = ltotold[chain];
	base_ltotnew = ltotnew[chain];

	base_ptotold = ptotold[chain];
	base_ptotnew = ptotnew[chain];
}



void screenOutput(int iter, int burnin, int eststeps, int burnin2)
{
	printf("\n ************ CHAIN %d ITERATION %d *****************", currentchain, iter);

	if (iter <= burnin)
	{
		printf("\n Phase: burn-in");
	}
	else if (iter < (burnin + eststeps))
	{
		printf("\n Phase: sampling");
	}
	else if (iter < (burnin + eststeps + burnin2))
	{
		printf("\n Phase: MLE search");
	}
	else
	{
		printf("\n Phase: MLE profile");
	}

	printf("\n likelihood_old \t %lf \t likelihood_new %lf \n", (double)ltotold[currentchain], (double)ltotnew[currentchain]);
	//printf("\n number OK %d \t",metr_number_ok);
	printf("\n -------------------------------------------- ");
	printf("\n acceptance ratio %lf \n samplesize %ld \n", (double)pbacc[currentchain] / (double)pbcount[currentchain], metr_number_ok[currentchain]);

	pbacc[currentchain] = 0.0;
	pbcount[currentchain] = 0.0;

	printf("\n --------------------------------------------");

	/* output 'display' parameters */
	printf("\n current parameter values:");
	printf("\n name\t\tcurrent_value\t(fixed)\t(num_accepted_changes, current_proposal_width)");
	param* params = current_para();
	for (int ii = 0; ii < paramcount; ii++)
	{
		if (params[ii].display > 0)
		{
			char* print_name = pprintname(ii);
			printf("\n %-14s %13lg\t(%d)\t(%d, %lf)",
				print_name,
				params[ii].pold,
				params[ii].fixed,
				params[ii].altt,
				params[ii].delta);
			delete[] print_name;
		}
	}

	printf("\n ********************************************* \n");
}

/*
 * Executes the likelihood functions. If the execution mode is set to parallel
 * the pfn_parallel_likelihood callback is excuted.
 */
void executeLikelihood()
{
#if DEBUG
	volatile DWORD dwStart = GetTickCount();
#endif

	if (parallel_mode == true)
	{
		// first set ltotnew to 0
		if (currentchain < 0)
		{
			base_ltotnew = 0.0;
		}
		else
		{
			ltotnew[currentchain] = 0.0;
			metr_number_ok[currentchain] = 0;
		}

		// now do calcs in parallel
#ifdef _OPENMP
		int proc = omp_get_max_threads();
#else
		int proc = 1;
#endif
		// work out sectioning
		int unit = likelihood_ubound - likelihood_lbound >= 0 ? 1 : -1;
		int scount = __max(1, __min(proc, unit*(likelihood_ubound - likelihood_lbound)));
		// scount==1 || scount <= abs(ubound-lbound)
		double delta = (double)(likelihood_ubound - likelihood_lbound) / (double)scount;
		// abs(delta) >= 1

		std::vector<double> tltot(scount);
		std::vector<long> tnumok(scount);

#pragma omp parallel for 
		for (int section = 0; section < scount; section++)
		{
			int start = likelihood_lbound + (int)(section*delta);
			int stop = likelihood_lbound + (int)((section + 1)*delta) - unit;
			if (section + 1 >= scount || unit*(stop - likelihood_ubound) > 0) stop = likelihood_ubound;
			double ltot = 0.0;
			long int numok = 0;

			pfn_parallel_likelihood(start, stop, &ltot, &numok); // this gets the plikelihood function to write into ltot and numok
			tltot[section] = ltot;
			tnumok[section] = numok;

		} // parallel for loop end

		if (currentchain < 0)
		{
			for (int section = 0; section < scount; section++)
			{
				base_ltotnew += tltot[section];
			}
		}
		else
		{
			for (int section = 0; section < scount; section++)
			{
				ltotnew[currentchain] += tltot[section];
				metr_number_ok[currentchain] += tnumok[section];
			}
		}
	}
	else
	{
		pfn_likelihood();
	}

#if DEBUG
	volatile DWORD dwStop = GetTickCount();
	printf("executeLikelihood() ticks: %d\n", dwStop - dwStart);
#endif
}

/*
* Initializes the likelihood for each chain.
*/
void init_likelihood()
{
	for (int c = 0; c < chaincount; ++c)
	{
		currentchain = c;
		executeLikelihood();
		lnlike_priors();
		ltotold[c] = ltotnew[c];
		ptotold[c] = ptotnew[c];
	}

	currentchain = -1;
}

/*
 * Initializes the accpetnace flags for each chain.
 */
void init_acceptance()
{
	delete[] accept_forced;
	accept_forced = new int[chaincount];
	CHECK(accept_forced != NULL, "Memory allocation failed.");

	delete[] reject_forced;
	reject_forced = new int[chaincount];
	CHECK(reject_forced != NULL, "Memory allocation failed.");

	for (int c = 0; c < chaincount; ++c)
	{
		accept_forced[c] = 0;
		reject_forced[c] = 0;
	}

}

/***************************************************************************/





/*
* Simple quicksort to be used in the calc_posterior function.
* We use a quick sort as it has a average runtime of O(n log n).
* The worst case of O(n^2) should never happen as the data is randomized allways.
*/
int divide(double data[], int left, int right)
{
	double x = data[right];

	int j = left - 1;

	for (int i = left; i < right; i++)
	{
		if (x < data[i])
		{
			j++;
			double tmp = data[j];
			data[j] = data[i];
			data[i] = tmp;
		}
	}

	data[right] = data[j + 1];
	data[j + 1] = x;

	return j + 1;
}

void quicksort(double data[], int left, int right)
{
	if (left < right)
	{
		int div = divide(data, left, right);
		quicksort(data, left, div - 1);
		quicksort(data, div + 1, right);
	}
}


double credible_endpoint(double level, std::vector<double>& sorted_samples, double limit_value, bool lower)
{
	// level - a required credibility from 0 to 1. E.g. 0.95 for 95% interval
	// sorted_samples[size] - an array of parameter samples sorted in descending order (larger to smaller)
	// limit_value - a maximum or minimum value of the parameter.
	// lower - do we need an upper or a lower endpoint

	size_t size = sorted_samples.size();
	assert(level >= 0 && level < 1);
	assert(size > 1);
	if (lower) { assert(limit_value < sorted_samples[size - 1]); }
	else { assert(limit_value > sorted_samples[0]); }

	double tail_level = 0.5*(1.0 - level); // tail_level < 0.5
	double tail_index = tail_level*(size + 1); // tail_index > 0 && tail_index < 0.5*size + 0.5, 0.5<0.5*size => tail_index<size
	size_t index = (size_t)tail_index;
	assert(index >= 0 && index < size);
	double h_index = tail_index - index;
	double p_index, p_index_1;
	if (lower)
	{
		if (index > 0)
		{
			p_index = sorted_samples[size - index];
			p_index_1 = sorted_samples[size - index - 1];
		}
		else // h_index==0
		{
			p_index = limit_value;
			p_index_1 = sorted_samples[size - 1];
		}
	}
	else // upper
	{
		if (index > 0)
		{
			p_index = sorted_samples[index - 1];
			p_index_1 = sorted_samples[index];
		}
		else // h_index==0
		{
			p_index = limit_value;
			p_index_1 = sorted_samples[0];
		}
	}
	return p_index*(1 - h_index) + p_index_1*h_index;
}


void calc_posteriors()
{
	int vlength = bayes_chain_length[currentchain]; // this is correct, position left from fill_bayestable will equal number of row in table
	param* params = current_para();
	if (vlength > 1)
	{
		std::vector<double> samples(vlength);
		for (int i = 3; i < paramcount + 3; i++)
		{
			// copy parameter value
			for (int j = 0; j < vlength; j++) samples[j] = bayestable[currentchain][j][i];
			std::sort(samples.begin(), samples.end(), std::greater<double>());
			double sum = 0;
			std::for_each(samples.begin(), samples.end(), [&sum](double v) {sum += v; });
			params[i - 3].bayes_mean = sum / vlength;
			params[i - 3].bcred_u95 = credible_endpoint(0.95, samples, params[i - 3].ub, false);
			params[i - 3].bcred_u68 = credible_endpoint(0.68, samples, params[i - 3].ub, false);
			params[i - 3].bcred_median = credible_endpoint(0.0, samples, params[i - 3].ub, false);
			params[i - 3].bcred_l68 = credible_endpoint(0.68, samples, params[i - 3].lb, true);
			params[i - 3].bcred_l95 = credible_endpoint(0.95, samples, params[i - 3].lb, true);
		}
	}
	else
	{
		for (int i = 0; i < paramcount; i++)
		{
			params[i].bcred_median = params[i].bayes_mean = 0.5 * (params[i].lb + params[i].ub);
			params[i].bcred_u95 = params[i].bcred_u68 = params[i].ub;
			params[i].bcred_l68 = params[i].bcred_l95 = params[i].lb;
		}
	}
}

/******************************************************/
void calc_MLE_intervals()
{
	double *mimin1 = new double[paramcount];
	double *mimin2 = new double[paramcount];
	double *mimax1 = new double[paramcount];
	double *mimax2 = new double[paramcount];
	int vlength;

	vlength = bayespos[currentchain];

	param* p = current_para();
	/* find MLE intervals (68% and 95%) */
	for (int ii = 0; ii < paramcount; ii++)
	{
		mimin1[ii] =
			mimin2[ii] =
			mimax1[ii] =
			mimax2[ii] = p[ii].MLE;
	}

	for (int pvcnt = 0; pvcnt < vlength; pvcnt++)
	{
		for (int ii = 0; ii < paramcount; ii++) {
			double tln = bayestable[currentchain][pvcnt][2];
			double ttpv = bayestable[currentchain][pvcnt][ii + 3];

			/* 68% */
			if ((ltotmax[currentchain] - tln) < 0.49447324) // qchisq(0.68,1)/2
			{
				// lowest value in this group
				if (ttpv < mimin1[ii])
					mimin1[ii] = ttpv;
				// highest value in this group
				if (ttpv > mimax1[ii])
					mimax1[ii] = ttpv;
			}

			/* 95% */
			if ((ltotmax[currentchain] - tln) < 1.92072941) // qchisq(0.95,1)/2
			{
				// lowest value in this group
				if (ttpv < mimin2[ii])
					mimin2[ii] = ttpv;
				// highest value in this group
				if (ttpv > mimax2[ii])
					mimax2[ii] = ttpv;
			}
		} // ii loop
	}

	param* params = current_para();
	for (int ii = 0; ii < paramcount; ii++)
	{
		params[ii].MLEl68 = mimin1[ii];
		params[ii].MLEl95 = mimin2[ii];
		params[ii].MLEu68 = mimax1[ii];
		params[ii].MLEu95 = mimax2[ii];
	}

	delete[] mimin1;
	delete[] mimin2;
	delete[] mimax1;
	delete[] mimax2;
}


/*
 * Prints general info, e.g. lotmax, aic, free params
 * details: (0) -- minimal details, (1) - full details
 */
void fprint_generalinfos(FILE *file, int chain, int details)
{

	const char* short_format = "metr_tag\t%s\n"
		"params\t%d\n"
		"free_params\t%d\n"
		"total_iterations\t%d\n"
		"sample_size\t%d\n"
		"max_likelihood\t%lf\n"
		"AIC\t%lf\n"
		"BIC\t%lf\n"
		"DIC\t%lf\n";

	const char* full_format = "metr_tag\t%s\n"
		"params\t%d\t\tTotal number of model parameters.\n"
		"free_params\t%d\t\tNumber of free (non-fixed) parameters.\n"
		"total_iterations\t%d\t\tNumber of itrations in burn-in and main sampling phases.\n"
		"sample_size\t%d\t\tUser-defined, usually equals to number of terms in likelihood.\n"
		"max_likelihood\t%lf\t\tMaximum log-likelihood reached.\n"
		"AIC\t%lf\t\tAkaike Information Criterion.\n"
		"BIC\t%lf\t\tSchwarz Information Criterion.\n"
		"DIC\t%lf\t\tDeviance Information Criterion.\n";

	const char* format = details > 0 ? full_format : short_format;

	fprintf(file, format,
		metr_tag,
		paramcount,
		chain < 0 ? base_freeparams : freeparams[chain],
		burnin + eststeps,
		chain < 0 ? metr_number_ok[bestchainID] : metr_number_ok[chain],
		chain < 0 ? base_ltotmax : ltotmax[chain],
		chain < 0 ? base_aic : aic[chain],
		chain < 0 ? base_bic : bic[chain],
		chain < 0 ? baseDic : dic[chain]);
}

/*
 * Prints the header or footer for the parameter table.
 * printFooter: (0) - no, print header; (1) - yes, print footer.
 */
void fprint_header(FILE *file, int printFooter)
{
	const char *header[] = {
		"ID","Numeric identifier",
		"name","text identifier together with numeric index for vector parameters if applicable.",
		"fixed?","(0) - free parameter, (1) - fixed parameter, (-1) - free parameter, but do not randomize it at the start.",
		"type", "(0) - real valued parameter, (1) - parameter that cannot be negative",
		"lb", "user-defined lower bound for parameter values",
		"ub", "user-difined upper bound for parameter values",
		"priorl68","user-defined lower bound of 68% confidence interval of prior parameter distribution; undefined if equals to -999",
		"priormode","most likely value (a.k.a. mode) of this parameter according to the prior; undefined if equals to -999",
		"prioru68","user-defined lower bound of 68% confidence interval of prior parameter distribution; undefined if equals to -999",
		"rootRhat", "estimated potential scale reduction across chains (see see Gilks et al., 'Markov Chain Monte Carlo in Practice', 1996, page 137)",
		"postl95", "lower bound of 95% confidence interval of posterior parameter distribution",
		"postl68", "lower bound of 68% confidence interval of posterior parameter distribution",
		"postmean", "mean value of posterior parameter distribution",
		"postu68", "upper bound of 68% confidence interval of posterior parameter distribution",
		"postu95", "upper bound of 95% confidence interval of posterior parameter distribution",
		"LP_l95", "lower bound of 95% band from likelihood profile method",
		"LP_l68", "lower bound of 68% band from likelihood profile method",
		"LP_MLE", "maximum likelihood estimate",
		"LP_u68", "upper bound of 68% band from likelihood profile method",
		"LP_u95", "upper bound of 95% band from likelihood profile method",
	};
	const int len = sizeof(header) / sizeof(header[0]);
	if (printFooter == 0)
	{
		for (int i = 0; i < len; i += 2)
		{
			if (i > 0) fputc('\t', file);
			fputs(header[i], file);
		}
		fputc('\n', file);
	}
	else
	{
		for (int i = 0; i < len; i += 2)
		{
			fprintf(file, "%s:\t%s\n", header[i], header[i + 1]);
		}
	}
}

/*
 * Prints the paramter specified by the given
 * index of the chain provided into the given file.
 */
void fprint_parameter(FILE *file, param *chain, int index)
{

	char* print_name = pprintname(index);
	fprintf(file, "%d\t%s\t", index, print_name);
	delete[] print_name;
	fprintf(file, "%d\t%d\t%lf\t%lf\t", chain[index].fixed, chain[index].type, chain[index].lb, chain[index].ub);
	if (chain[index].prioryn > 0 && chain[index].type == 0)
	{
		fprintf(file, "%lf\t%lf\t%lf\t", chain[index].priormean - chain[index].priorsdev, chain[index].priormean, chain[index].priormean + chain[index].priorsdev);
	}
	if (chain[index].prioryn > 0 && chain[index].type == 1)
	{
		fprintf(file, "%lf\t%lf\t%lf\t", exp(log(chain[index].priormean) - chain[index].priorsdev), chain[index].priormean, exp(log(chain[index].priormean) + chain[index].priorsdev));
	}
	if (chain[index].prioryn <= 0)
	{
		fprintf(file, "-999.0 \t-999.0 \t-999.0 \t");
	}
	// convergence stat
	fprintf(file, "%lf\t", chain[index].rootRhat);
	fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t", chain[index].bcred_l95, chain[index].bcred_l68, chain[index].bayes_mean, chain[index].bcred_u68, chain[index].bcred_u95);
	fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t", chain[index].MLEl95, chain[index].MLEl68, chain[index].MLE, chain[index].MLEu68, chain[index].MLEu95);
	fputc('\n', file);
}

/*
 * Prints <name>_final_out.txt and <name>_params_out.txt files
 */
void final_metro_output(FILE *mmfile, FILE *mmfile2)
{

	if (out_summary)
	{
		// foreach of the chains run calc_IC
		fprintf(mmfile,
			"========================================\n"
			"Selected Chain is Chain %d \n"
			"========================================\n\n"
			, bestchainID);

		fprint_generalinfos(mmfile, -1, 1);
		fputc('\n', mmfile);

		fprint_header(mmfile, 0);	// print header	
		for (int p = 0; p < paramcount; p++)
		{
			fprint_parameter(mmfile, current_para(-1), p); // print table
		}
		fputc('\n', mmfile);
		fprint_header(mmfile, 1);	// print footer
		fputc('\n', mmfile);

		if (chaincount > 0) fprintf(mmfile,
			"========================================\n"
			"Other Chains\n"
			"========================================\n");

		for (int c = 0; c < chaincount; c++) if (c != bestchainID)
		{
			fprintf(mmfile, "\nChain %i\n\n", c);

			fprint_generalinfos(mmfile, c, 0);
			fputc('\n', mmfile);

			fprint_header(mmfile, 0);
			for (int p = 0; p < paramcount; p++)
			{
				fprint_parameter(mmfile, current_para(c), p);
			}
		}
		fclose(mmfile);
	}

	// output for final params output file, same as best chain

	if (out_params)
	{
		fprint_header(mmfile2, 0);
		for (int p = 0; p < paramcount; p++)
		{
			fprint_parameter(mmfile2, current_para(-1), p);
		}
		fputc('\n', mmfile2);
		fclose(mmfile2);
	}
}

/******************************************************/
void params_draw_random_vector()
{
	int chain = currentchain;
	if (chain < 0) chain = bestchainID;
	int vlength = bayes_chain_length[chain];
	param* params = current_para();
	// choose random number down the vlength
	int tnum = ((int)(drand()*vlength));

	for (int mm = 0; mm < paramcount; mm++)
	{
		params[mm].value = bayestable[currentchain][tnum][mm + 3];
	}
}

/******************************************************/
// returns 0 if successful
int params_from_bayes_list(int chain, int index)
{
	if ((chain < 0) || (chain >= chaincount)) return 1;
	int vlength = bayespos[chain];

	if (vlength > bayes_chain_length[chain])
	{
		vlength = bayes_chain_length[chain]; // then going to loop j= 0 to vlength only
	}

	if ((index < 0) || (index >= vlength)) return 1;

	param* params = current_para();
	for (int mm = 0; mm < paramcount; mm++)
	{
		params[mm].value = bayestable[chain][index][mm + 3];
	}
	return 0;
}


/******************************************************/
int params_from_bayes_list(int index)
{
	if (currentchain < 0) return params_from_bayes_list(bestchainID, index);
	return params_from_bayes_list(currentchain, index);
}


/*******************************************************/
void params_set_to_posterior_mean()
{
	param* params = current_para();
	for (int ii = 0; ii < paramcount; ii++)
	{
		params[ii].value = params[ii].bayes_mean;
	}
}

/*
 * Calculates no. free params, and AIC, BIC, and DIC.
 */
void calc_IC()
{
	int *tempFreeP;
	double *tempAic;
	double *tempBic;
	double *tempDic;
	double *templtotmax;


	param* params = current_para();
	if (currentchain < 0)
	{
		tempFreeP = &base_freeparams;
		tempAic = &base_aic;
		tempBic = &base_bic;
		tempDic = &baseDic;
		templtotmax = &base_ltotmax;

		for (int i = 0; i < paramcount; ++i)
		{
			if (params[i].fixed < 1)
			{
				(*tempFreeP)++;
			}
		}
	}
	else
	{
		tempFreeP = &freeparams[currentchain];
		tempAic = &aic[currentchain];
		tempBic = &bic[currentchain];
		tempDic = &dic[currentchain];
		templtotmax = &ltotmax[currentchain];

		for (int i = 0; i < paramcount; ++i)
		{
			if (params[i].fixed < 1)
			{
				(*tempFreeP)++;
			}
		}
	}

	/* aic */
	(*tempAic) = ((-2.0) * (*templtotmax)) + 2.0 * ((double)*tempFreeP);
	/* bic */
	(*tempBic) = ((-2.0) * (*templtotmax)) + log((double)metr_number_ok[currentchain]) * ((double)*tempFreeP);

	/* dic */
	// DIC computation reworked
	CHECK(currentchain >= 0, "calc_IC cannot be called out of chain");
	double dmean_sum1 = 0.0;

	int n = bayes_chain_length[currentchain];
	for (int i = 0; i < n; ++i)
	{
		dmean_sum1 += (-2.0) * bayestable[currentchain][i][1]; /* not sure about this when using priors -- check */
	}
	dmean_sum1 /= n;

	params_set_to_posterior_mean();
	executeLikelihood();

	// lnlike_priors(); /* not sure about this -- check */

	double datmean = (-2.0) * (ltotnew[currentchain] + ptotnew[currentchain]); /* not sure about this when using priors -- check */
	// pD = dmean_sum1 - datmean;
	(*tempDic) = (-1.0) * datmean + (2.0) * dmean_sum1;

}

/******************************************************/
void calc_convergence()
{
	param* base_para = current_para(-1);
	if (chaincount > 1)
	{

		int m = chaincount;
		double* psi_i_mean = new double[m];
		double* si_squared = new double[m];

		convergence_RhatMean = 0.0;
		int Rhatmean_count = 0;

		int n = bayes_chain_length[0];
		for (int c = 1; c < chaincount; c++)
			if (n > bayes_chain_length[c])
				n = bayes_chain_length[c];

		// loop over params
		for (int mm = 0; mm < paramcount; mm++)
		{
			if (base_para[mm].fixed < 1)
			{
				// first calculate psi_i_mean and si_squared for this param for each chain i, also calculate grand mean psi_mean
				double psi_mean = 0.0;
				for (int i = 0; i < m; i++) // loop over chains
				{
					psi_i_mean[i] = 0.0;
					// loop over samples to calculate mean
					for (int j = 0; j < n; j++)
					{
						psi_i_mean[i] += bayestable[i][j][mm + 3];
					}
					psi_i_mean[i] /= n;

					// update grand mean i.e. mean over all chains
					psi_mean += psi_i_mean[i];

					// squared differences within this chain
					si_squared[i] = 0.0;
					for (int j = 0; j < n; j++)
					{
						si_squared[i] += pow(bayestable[i][j][mm + 3] - psi_i_mean[i], 2);
					}
					si_squared[i] /= n - 1;
				} // end of loop over chains

				psi_mean /= m;

				// calc capB for this param
				double capB = 0.0;
				for (int i = 0; i < m; i++)
				{
					capB += pow(psi_i_mean[i] - psi_mean, 2);
				}
				capB *= n / (m - 1.0);

				// capW
				double capW = 0.0;
				for (int i = 0; i < m; i++)
				{
					capW += si_squared[i];
				}
				capW /= m;

				// varhat
				double varhat = ((n - 1)*capW + capB) / n;

				// rootRhat
				double rootRhat = sqrt(varhat / capW);

				for (int c = -1; c < chaincount; c++)
					current_para(c)[mm].rootRhat = rootRhat;

				convergence_RhatMean += base_para[mm].rootRhat;
				Rhatmean_count += 1;
			} // if para not fixed
			else
			{
				base_para[mm].rootRhat = -999.0;
			}


		} // loop over params

		convergence_RhatMean /= Rhatmean_count;

		delete[] psi_i_mean;
		delete[] si_squared;

	}
	else
	{
		for (int mm = 0; mm < paramcount; mm++)
		{
			base_para[mm].rootRhat = -999.0;
		}
		convergence_RhatMean = -999.0;
	}

}

/******************************************************/
void force_accept()
{
	accept_forced[currentchain] = 1;
}

/*******************************************************/
void force_reject()
{
	reject_forced[currentchain] = 1;
}

/*
* TODO: Some meaningfull description here
*/
void params_set_to_old()
{
	param* p = current_para();
	for (int i = 0; i < paramcount; i++)
	{
		p[i].value = p[i].pold;
	}
}

/*******************************************************/
int inburnin()
{
	if (metr_k < burnin) return 1;
	return 0;

}

void runmcmc(int tburnin, int teststeps, int tburnin2, int tmleexp)
{

	//time_t start;
	//time_t end;
	//time_t startParamLoop;
	//time_t endParamLoop;

	//time(&start);
	//volatile DWORD dwStart = GetTickCount();


	FILE *mfile = NULL;
	FILE *MLEfile = NULL;
	FILE *bayesfile2 = NULL;
	FILE *outfile = NULL;
	FILE *outfile2 = NULL;

	CHECK(tburnin > 0, "First parameter to runmcmc MUST be greater than zero");
	CHECK(teststeps >= bayes_step, "Second parameter to runmcmc cannot be less than 100");
	/* set burnin etc. */
	burnin = tburnin;
	burnin2 = tburnin;
	eststeps = teststeps;
	mleexp = tmleexp;

	if (out_mle)
	{
		MLEfile = workspace_fopen(out_mle_name, "w");
		if (MLEfile != NULL) {
			// do header 
			// chain, iteration number
			fprintf(MLEfile, "chain\titeration\tposterior\tlikelihood");
			// parameter names
			for (int ii = 0; ii < paramcount; ii++)
			{
				char* print_name = pprintname(ii);
				fprintf(MLEfile, "\t%s", print_name);
				delete[] print_name;
			}
			fprintf(MLEfile, "\n");
		}
	}

	if (out_bayes)
	{
		bayesfile2 = workspace_fopen(out_bayes_name, "w");
		if (bayesfile2 != NULL) {
			// do header for bayes list files
			// chain, iteration number
			fprintf(bayesfile2, "chain\titeration\tposterior\tlikelihood");
			// parameter names
			for (int ii = 0; ii < paramcount; ii++)
			{
				char* print_name = pprintname(ii);
				fprintf(bayesfile2, "\t%s", print_name);
				delete[] print_name;
			}
			fprintf(bayesfile2, "\n");
		}
	}

	if (out_summary)
	{
		outfile = workspace_fopen(out_summary_name, "w");
	}

	if (out_params)
	{
		outfile2 = workspace_fopen(out_params_name, "w");
	}


	printf("\nbeginning metropolis \n");
	printf("number params \t %d \n", paramcount);

	// Initialization
	init_tots();
	init_chains(chaincount);
	init_acceptance();
	init_counters();
	init_bayestable(((burnin + teststeps) / bayes_step) - ((burnin + bayes_step - 1) / bayes_step) + 1, tburnin2 + tmleexp);
	// also keep track of length of the truly 'bayesian' part of this table
	clear_altlists();

	/* some little calcs needed for later on */
	delete[] numalt;
	numalt = countFreeParams();

	std::unique_ptr<double> talt(calcTAlt(numalt));

	initStepSizes();
	initRandomValues();
	init_likelihood();
	init_paramerniformation();
	init_freeparams();
	init_ic();

	if (out_chain)
	{
		mfile = workspace_fopen(out_chain_name, "w");
		// do header 
		// chain, iteration number
		fprintf(mfile, "chain\titeration\tposterior\tlikelihood");
		// parameter names
		for (int ii = 0; ii < paramcount; ii++)
		{
			char* print_name = pprintname(ii);
			fprintf(mfile, "\t%s", print_name);
			delete[] print_name;
		}
		fprintf(mfile, "\n");
		for (int chain = 0; chain < chaincount; chain++)
		{
			param* params = current_para(chain);
			fprintf(mfile, "%d\t%d\t%lf\t%lf", chain, metr_k, ltotnew[chain] + ptotnew[chain], ltotnew[chain]);
			for (int i = 0; i < paramcount; i++)
			{
				fprintf(mfile, "\t%lf", params[i].value);
			}
			fprintf(mfile, "\n");
		}
		fclose(mfile);
	}

	/***********************************************************************/
	/*****  metropolis loop  ***********************************************/
	/***********************************************************************/

	/* loop over chains */
#ifdef _OPENMP
	omp_set_nested(1);
#endif
#pragma omp parallel for
	for (int chain = 0; chain < chaincount; chain++)
	{
		// chain specific variables
		bool accept_flag = false; // flag at least one accept since last sample
		for (metr_k = 1; metr_k <= (tburnin + teststeps + tburnin2 + tmleexp); metr_k++)
		{

			currentchain = chain;
			param* chain_params = current_para(chain);
			/* --------- choose which parameters to change ----------- */
			/* select parameters to change */
			/* most of the time change only one, two or three params */

			if (drand() < 0.670)
			{
				for (int i = 0; i < paramcount; i++)
				{
					chain_params[i].alt = 0;
				}

				int rnd = 0;
				/* choose one param to alter */
				do //BUG: This might just run forever, needs to be refactored
				{
					rnd = (int)(drand() * paramcount);
					if (chain_params[rnd].fixed < 1 && chain_params[rnd].delay < metr_k)
					{
						chain_params[rnd].alt = 1;
						break;
					}
				} while (true);

				/* alter parameters close by? */
				if ((rnd - 1) >= 0 && drand() < 0.500)
				{
					if (chain_params[rnd - 1].fixed < 1 && chain_params[rnd - 1].delay < metr_k)
						chain_params[rnd - 1].alt = 1;
				}

				if ((rnd + 1) < paramcount && drand() < 0.500)
				{
					if (chain_params[rnd + 1].fixed < 1 && chain_params[rnd + 1].delay < metr_k)
						chain_params[rnd + 1].alt = 1;
				}
			}
			else
			{
				/* change many parameters at once */
				/* draw prob change for this iteration */
				double palt = talt.get()[chain] * exp(4.0*(drand() - 0.50));

				if (palt < (0.10*talt.get()[chain]))
					palt = 0.10*talt.get()[chain];

				if (palt > 0.99)
					palt = 0.99;

				int calt = 0;
				do //BUG: Potential infinite runtime, needs to be addressed
				{
					for (int ii = 0; ii < paramcount; ii++)
					{
						if (chain_params[ii].fixed < 1 && chain_params[ii].delay < metr_k && drand() < palt)
						{
							chain_params[ii].alt = 1;
							calt++;
						}
						else
						{
							chain_params[ii].alt = 0;
						}
					}
				} while (calt == 0);
			} // if changing one-to-many parameters

			/* ---------- change parameter values, i.e. make the 'jump' ----------------- */

			for (int ii = 0; ii < paramcount; ii++)
			{
				if (chain_params[ii].alt > 0)
				{
					do //BUG: This runs potentially for ever, nneds o be addressed
					{
						// mpara[ii].padd = ((drand48()-0.50)*2.0)*mpara[ii].delta;
						chain_params[ii].padd = normal_draw(0.0, chain_params[ii].delta);

						if (chain_params[ii].type > 0)
							chain_params[ii].value = chain_params[ii].pold*exp((double)chain_params[ii].padd);
						else
							chain_params[ii].value = chain_params[ii].pold + chain_params[ii].padd;

						if (chain_params[ii].value<chain_params[ii].ub && chain_params[ii].value>chain_params[ii].lb)
							break;
						else
							chain_params[ii].value = chain_params[ii].pold;

					} while (true);

				} // if altered

			} // ichg loop

			/* ---------- MLE tails exploration ------------- */
			/* start of this phase */
			if (metr_k == ((burnin + eststeps) + 1))
			{
				ASSERT(bayespos[chain] <= bayestable_length, "bayestable_length improperly computed")
					bayes_chain_length[chain] = bayespos[chain];
				for (int ii = 0; ii < paramcount; ii++)
				{
					if (chain_params[ii].fixed < 1)
					{
						chain_params[ii].value = chain_params[ii].MLE;
						// mpara[ii].delta *= 0.10; 
					}
				}
			}


#if 1
			/* set back to MLE values (i.e. in second burning period?) */
			if (metr_k > (burnin + eststeps) && metr_k < (burnin + eststeps + burnin2) && metr_k % 50 == 0)
			{
				for (int ii = 0; ii < paramcount; ii++)
				{
					if (chain_params[ii].fixed < 1)
					{
						/* i.e. keeps resetting to MLE */
						chain_params[ii].value = chain_params[ii].MLE;
					}
				}
			}
#endif


			/* ------------ calc new lnlike ------------- */
			executeLikelihood();
			lnlike_priors();

			/* ---------- compare "new" likelihood to "max" likelihood to see if we have reached the best fit yet (or not) ------------ */
			if (ltotnew[chain] > ltotmax[chain] || metr_k < 2)
			{
				ltotmax[chain] = ltotnew[chain];

				for (int ii = 0; ii < paramcount; ii++)
				{
					chain_params[ii].MLE = chain_params[ii].value;
				}

				/* update altt for each para */
				if (metr_k >= 2 && reject_forced[chain] != 1)
				{
					for (int ii = 0; ii < paramcount; ii++)
						if (chain_params[ii].alt>0)
							chain_params[ii].altt++;
				}

			} // if ltotnew>ltotmax

			/* ------------------ screen output ------------------- */

			if ((console_out_step > 0) && (metr_k%console_out_step == 0))
			{
#pragma omp critical
			{
				screenOutput(metr_k, burnin, eststeps, burnin2);
			}
			}

			/* ------------- output to MLE list file for use later on? ----------------------------- */
			/* ------------- nb only do this when past end of sampling phase ----------------------- */
			if ((metr_k > (tburnin + teststeps)) && (ltotmax[chain] - ltotnew[chain]) <= 2.00)
			{
				// output to table held in memory (actually just add to bayestable)
				fill_bayestable(chain, metr_k);
			}

			/* --------- compare new to old and accept or reject -- METROPOLIS CRITERION IS IN HERE ------------- */
			bool accept = false;

			if (((ltotnew[chain] + ptotnew[chain]) > (ltotold[chain] + ptotold[chain])) || accept_forced[currentchain] > 0)
			{
				/* always accept if ltotnew is higher, or we have been forced to accept */
				accept = true;
			}
			else
			{
				/* otherwise accept with probablity based on difference between ltotnew and ltotold */
				double rndnum = (double)drand();
				if (rndnum < 0.00000010)
				{
					rndnum = 0.00000010;
				}
				if (rndnum > 0.99999990)
				{
					rndnum = 0.99999990;
				}

				double dlik = (double)((ltotnew[chain] + ptotnew[chain]) - (ltotold[chain] + ptotold[chain]));
				if (log(rndnum) < dlik)
				{
					accept = true;
				}
				else
				{
					accept = false;
				}
			} // end of acceptance criterion bit

			// did we force reject?
			if (reject_forced[currentchain] > 0)
			{
				accept = false;
			}

			/* for keeping track of overall acceptance ratio */
			if (accept_forced[chain] != 1 && reject_forced[chain] != 1)
			{
				pbcount[chain] += 1.0;
				if (accept)
				{
					pbacc[chain] += 1.0;
				}

			}
			/* --------- act on acceptance --------- */
			if (accept)
			{
				/* update likelihoods memory */
				ltotold[chain] = ltotnew[chain];
				ptotold[chain] = ptotnew[chain];

				/* update "old" parameters */
				for (int ii = 0; ii < paramcount; ii++) {

					if (chain_params[ii].fixed < 1)
						chain_params[ii].pold = chain_params[ii].value;

					if (chain_params[ii].alt > 0)
					{
						/* running totals */
						if (accept_forced[chain] < 1)
						{
							chain_params[ii].runalt++;
							chain_params[ii].runacc++;
						}

						/* set decision array */
						chain_params[ii].decis = 2;
					} // if alt

				} // ii loop
				accept_flag = true;
			}
			else
			{
				/* --------------- act on rejection --------------- */

				/* update likelihood memory */
				ltotnew[chain] = ltotold[chain];
				ptotnew[chain] = ptotold[chain];

				for (int ii = 0; ii < paramcount; ii++)
				{
					/*return to "old" parameters*/
					if (chain_params[ii].fixed < 1)
						chain_params[ii].value = chain_params[ii].pold;

					if (chain_params[ii].alt > 0) {

						/* running totals */
						if (reject_forced[chain] < 1)
						{
							chain_params[ii].runalt++;
						}

						/* set decision array */
						chain_params[ii].decis = 1;

					} // if param altered
				} // ii loop
			} // if reject

			/* adjust jump steps? only if in burnin */
			if (metr_k < burnin)
			{
				for (int ii = 0; ii < paramcount; ii++)
				{
					if (chain_params[ii].runalt == 20 && chain_params[ii].fixed < 1)
					{
						chain_params[ii].runalt = 0;

						if (chain_params[ii].runacc < 5)
						{
							/* decrease temperature */
							chain_params[ii].delta *= 0.80;
						}
						if (chain_params[ii].runacc > 5)
						{
							/* increase temperature */
							chain_params[ii].delta *= 1.20;
						}
						chain_params[ii].runacc = 0;

						if (chain_params[ii].type > 0) {
							if (chain_params[ii].delta < 0.010)
								chain_params[ii].delta = 0.010;
							if (chain_params[ii].delta > 10.0)
								chain_params[ii].delta = 10.0;
						}
						else {
							if (chain_params[ii].delta < (0.0010*(chain_params[ii].ub - chain_params[ii].lb)))
								chain_params[ii].delta = (0.0010*(chain_params[ii].ub - chain_params[ii].lb));
							if (chain_params[ii].delta > (0.50*(chain_params[ii].ub - chain_params[ii].lb)))
								chain_params[ii].delta = (0.50*(chain_params[ii].ub - chain_params[ii].lb));
						}
					} // if runalt == 20
				} // ii loop
			} // if still in first burnin

			/* keeping track of decisions thing */
			for (int ii = 0; ii < paramcount; ii++)
				if (chain_params[ii].alt < 1)
					chain_params[ii].decis = 0;

			if (metr_k % bayes_step == 0 && metr_k >= tburnin  && metr_k <= (tburnin + teststeps) && accept_flag)
			{
				// write sample into memory
				fill_bayestable(chain, metr_k);
				accept_flag = false;
			} // if outputting sample

			/* -------------- write to another file that logs everything (this file not used any more) ------------- */
			if (accept && out_chain)
			{
#pragma omp critical
			{
				mfile = workspace_fopen(out_chain_name, "a");
				// output chain, iteration
				fprintf(mfile, "%d\t%d\t%lf\t%lf", chain, metr_k, ltotnew[chain] + ptotnew[chain], ltotnew[chain]);
				// parameter values
				param* p = current_para();
				for (int i = 0; i < paramcount; i++)
				{
					fprintf(mfile, "\t%lf", p[i].value);
				}
				fprintf(mfile, "\n");
				fclose(mfile);
			}

			}

			/* unforce any forced accepts, rejects */
			accept_forced[currentchain] = 0;
			reject_forced[currentchain] = 0;

		} // metr_k
	} // chain


	printf("\n ****************************************************************");
	printf("\n metropolis loop finished, beginning final calcs and output");
	printf("\n ****************************************************************\n");

	if (bayesfile2 != NULL)
	{
		for (int cc = 0; cc < chaincount; cc++)
		{
			// output bayes list allchains file
			for (int bb = 0; bb < bayes_chain_length[cc]; bb++)
			{
				// output chain, iteration
				fprintf(bayesfile2, "%d\t%d\t%lf\t%lf", cc, (int)bayestable[cc][bb][0], bayestable[cc][bb][1], bayestable[cc][bb][2]);
				// parameter values 
				for (int i = 0; i < paramcount; i++)
				{
					fprintf(bayesfile2, "\t%lf", bayestable[cc][bb][i + 3]);
				}
				fprintf(bayesfile2, "\n");
			}
		}
		fclose(bayesfile2);
	}
	if (MLEfile != NULL)
	{
		// output MLE list allchains file
		for (int cc = 0; cc < chaincount; cc++)
		{
			for (int bb = bayes_chain_length[cc]; bb < bayespos[cc]; bb++)
			{
				// output chain, iteration
				fprintf(MLEfile, "%d\t%d\t%lf\t%lf", cc, (int)bayestable[cc][bb][0], bayestable[cc][bb][1], bayestable[cc][bb][2]);
				// parameter values 
				for (int i = 0; i < paramcount; i++)
				{
					fprintf(MLEfile, "\t%lf", bayestable[cc][bb][i + 3]);
				}
				fprintf(MLEfile, "\n");
			}
		}
		fclose(MLEfile);
	}

	calc_convergence();

	currentchain = -1;

#pragma omp parallel for
	for (int chain = 0; chain < chaincount; chain++)
	{
		currentchain = chain;


		/* calc posteriors from samples */

		calc_posteriors();

		/* calc MLE intervals */

		calc_MLE_intervals();

		/* set params to bayesian averages */
		param* params = current_para(chain);
		for (int ii = 0; ii < paramcount; ii++)
		{
			params[ii].value = params[ii].bayes_mean;
		}

		/* final calc of lnlike to set things up nicely */
		for (int ii = 0; ii < paramcount; ii++)
		{
			params[ii].alt = 1;
		}

		executeLikelihood();
		lnlike_priors();

	} // chain
	currentchain = -1;

#pragma omp parallel for
	for (int chain = 0; chain < chaincount; chain++)
	{
		currentchain = chain;
		calc_IC();
	}
	currentchain = -1;

	consolidateChains();

	/* final output */
	final_metro_output(outfile, outfile2);

	printf("\n ****************************************************************");
	printf("\n metropolis all done");
	printf("\n no. chains: %d, runmcmc phases: %d,%d,%d,%d\n overall convergence (close to 1.0 is good, >1.2 is worse): %lf",
		chaincount, burnin, teststeps, tburnin2, tmleexp, convergence_RhatMean);
	printf("\n ****************************************************************\n");

}
