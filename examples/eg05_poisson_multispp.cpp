#include "preamble.h"
#if MODEL == 5

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "filzbach.h"

void pause(){PAUSE}

/************************************************************/
/* function headers                                         */
/************************************************************/

void fake_data();
void read_data();
void setup_parameters();
void fit_model();
void final_output();

/************************************************************/

void fake_data()
{
	/* use this to create fake data to test your analysis, before you try real data */

	/* how many data? */
	int numdata = 999;

	/* how many species? */
	int numsp = 10;
	/* generate random means. */
	/* do we get these parameters back out again? */
	table_create("trueparams");
	table_addcolumn("trueparams","spid");
	table_addcolumn("trueparams","lambda");
	for(int ss = 0; ss < numsp; ss++)
	{
		table_writevalue("trueparams", "spid", ss, ss);
		table_writevalue("trueparams", "lambda", ss, 1.0+random(0.0,9.0));
	}

	table_create("fakedata");
	table_addcolumn("fakedata","spid");
	table_addcolumn("fakedata","count");
	for(int ii = 0; ii < numdata; ii++)
	{
		/* generate random species */
		int spp = random_integer(0,numsp-1);
		/* draw random value for y from poisson distribution with appropriate mean and sdev */
		double true_lambda = table_getvalue("trueparams", "lambda", spp);
		int y = poisson_draw(true_lambda);
		/* write to fake data table */
		table_writevalue("fakedata", "spid", ii, spp);
		table_writevalue("fakedata", "count", ii, y);
	}

	table_output("fakedata","./workspace/eg05_poisson_multispp_fakedata.txt");
	table_output("trueparams","./workspace/eg05_poisson_multispp_trueparams.txt");

	return;
}

/************************************************************/
void read_data()
{
	table_read("mydata","./workspace/eg05_poisson_multispp_fakedata.txt",2);

	PAUSE
	return;
}

/************************************************************/

void setup_parameters()
{
	/* first, figure out what is the number of species */
	int numdata = table_numrows("mydata");
	int numsp = (int)table_getcolumnmax("mydata","spid") + 1;

	/* each line defines one new parameter for our model */
	parameter_create_vector("lambda", 0.0010, 200.0, 100.0, 1, 0, 1, numsp);

	/* set parameters for the hierarchical distribution? */
	parameter_create("lambda_mean",0.0010,200.0,100.0,1,0,1);
	parameter_create("lambda_sdev",0.0010,2.0,0.20,1,0,1);
	
	parameter_showall();

	PAUSE
	return;
}  

/************************************************************/

void likelihood()
{
	/* writes the log-likelihood given current parameter values */

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* get model parameters from list held by metropolis header */

	/* loop over data */
	int numdata = table_numrows("mydata");
	int numsp = 0;
	for(int ii = 0; ii < numdata; ii++)
	{
		/* get observed count and species id */
		int count = (int)table_getvalue("mydata","count",ii);
		int spp = (int)table_getvalue("mydata","spid",ii);
		if (spp > numsp) numsp = spp;
		/* get parameter for this species */
		double prob1 = poisson_density(count, cv("lambda",spp));
		inc_metr_ltotnew(log(prob1));
		inc_metr_number_ok(1);
	}
	numsp += 1;

	/* loop over parameters for hierarhical modelling */
	double grandmean = cv("lambda_mean");
	double grandsdev = cv("lambda_sdev");
	for(int spp = 0; spp < numsp; spp++)
	{
		/* get param value for this species */
		double mean = cv("lambda",spp);
		/* calc prob1 */
		double prob1 = normal_density(log(mean),log(grandmean),grandsdev);
		inc_metr_ltotnew(log(prob1));
	}
	
	return;
}

/************************************************/

void final_output()
{
	/* create link to file */
	char fname[100];

	/* create file name for output -- outp is the 'path' */
	get_filzbach_path(fname, 100);
	/* now your bit to add to the path */
	strcat(fname,"_my_output.txt");

	/* print internal table out to file with name taken from above */
	table_output("mydata",fname);

	return;
}

/* ************************************************* */
/* The control function that calls everything else.  */
/* ************************************************* */

int main()
{
	atexit(pause);
  // set the likelihood function pointer
  pfn_likelihood = &likelihood;

  initialize_filzbach();
  name_analysis("eg05_poisson_multispp_00");
  fake_data(); // either fake data (to test routine) or read data (to do for real)
  read_data();
  setup_parameters();
  //set_chains(3);
  runmcmc(5000, 5000, 5000, 5000);
  final_output();
}

#endif

