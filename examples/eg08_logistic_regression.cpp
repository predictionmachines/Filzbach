#include "preamble.h"
#if MODEL == 8

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
	/* do we get these parameters back out again? */
	double true_k0 = -2.10;
	double true_alpha = 0.025;

	/* how many data in all? */
	int numdata = 1000;

	/* create internal table to hold data, and model predictions (written in later) */
	table_create("fakedata");
	table_addcolumn("fakedata","x");
	table_addcolumn("fakedata","y");

	for(int ii = 0; ii < numdata; ii++)
	{
		/* draw random value for x from uniform distribution */
		double x = random(0.0,1.0)*100.0;
		/* calculte the 'k' value for this x */
		double k = true_k0 + true_alpha*x;
		/* convert to prob of 1 */
		double prob = logistic(k);
		// same as prob = 1.0 / (1.0 + exp((-1.0)*k)));
		/* draw 1 or 0 at random, with probability of 1 being 'prob' */
		double y;
		if (random(0.0,1.0)<=prob)
			y = 1.0;
		else
			y = 0.0;
		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata","x",ii,x);
		table_writevalue("fakedata","y",ii,y);		
	}

	table_output("fakedata","./workspace/eg10_logisticregression_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read table from external file into internal file called 'mydata' */
	table_read("mydata","./workspace/eg10_logisticregression_fakedata.txt",2);

	/* in this case, we are going to write into our internal table, the model predictions */
	/* so we need an extra columnd for that */
	table_addcolumn("mydata","predprob");

	PAUSE
		return;

}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("k0",-5.0,5.0,0.0,0,0,1);
	parameter_create("alpha",-0.20,0.20,0.0,0,0,1);
	parameter_create("beta",-0.050,0.050,0.0,0,0,1);

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

	/* loop over data */
	int numdata = table_numrows("mydata");
	for(int ii = 0; ii < numdata; ii++)
	{
		double x = table_getvalue("mydata","x",ii);
		double y = table_getvalue("mydata","y",ii);
		/* calc k */
		double k = cv("k0") + cv("alpha")*x + cv("beta")*x*x; 
		// not model we are fitting contains term that was not present in 'reality', 
		// i.e. the fake data

		/* prob of 1 */
		double predprob = logistic(k);
		/* compare to real data */
		double prob1;
		if(y>0.0)
			prob1 = predprob;
		else
			prob1 = 1.0 - predprob;

		if(prob1<0.000000010)
			prob1 = 0.000000010;

		// write prediction into data structure 
		// (nb would not be allowable under OMP parallel chains: serial chains only)
		table_writevalue("mydata","predprob",ii,predprob);

		// add log prob to sum
		inc_metr_ltotnew(log(prob1));
		inc_metr_number_ok(1);
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
	name_analysis("eg10_logisticregression_00");
	fake_data();
	read_data();
	setup_parameters();
	runmcmc(5000, 50000, 5000, 5000);
	final_output();
}

#endif

