#include "preamble.h"
#if MODEL == 12

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "filzbach.h"

void pause(){PAUSE}

/************************************************************/
/* global variables                                         */
/************************************************************/

/* some things needed for book-keeping (see below for how each one is used) */
int numdata;

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
	double true_m0 = 2.0, true_beta=0.010, true_gph=0.40, true_gdepth=0.00, true_gmoss=0.20; // nb -- no effect of depth
	double phref=7.0, depthref=10.0, mossref=0.0;
	double ph, depth, moss, lambda, count;
	int ii;

	/* create internal table to hold data, and model predictions (written in later) */
	table_create("fakedata");
	table_addcolumn("fakedata","ph");
	table_addcolumn("fakedata","depth");
	table_addcolumn("fakedata","moss");
    table_addcolumn("fakedata","count");

	/* how many data in all? */
	numdata = 999;

	for(ii=0 ; ii<numdata ; ii++)
	{
		/* draw fakedata values for factors: ph, depth, moss */
		depth = random(0.0,20.0);
		ph = 5.0 + 2.0*(depth/20.0) + random(0.0,2.0);
		
		if(random(0.0,1.0)<0.50)
			moss = 1.0;
		else
			moss = 0.0;

		/* calculte lambda (=average seedling count) for these conditions */
		lambda = true_beta + true_m0*exp( true_gph*(ph-phref) + true_gdepth*(depth-depthref) + true_gmoss*(moss-mossref));
		/* draw fakedata no. of seedlings from poisson with this lambda */
		count = poisson_draw(lambda);

		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata","moss",ii,moss);
		table_writevalue("fakedata","ph",ii,ph);
		table_writevalue("fakedata","depth",ii,depth);	
		table_writevalue("fakedata","count",ii,count);		
	}

	table_output("fakedata","./workspace/EgDataz2_seedrain.txt");

	return;
}

/************************************************************/
void read_data()
{
    /* read table from external file into internal file called 'mydata' */
    table_read("mydata","./workspace/EgDataz2_seedrain.txt",4);
	numdata = table_numrows("mydata");

	/* in this case, we are going to write into our internal table, the model predictions */
	/* so we need an extra columnd for that */
	table_addcolumn("mydata","predlambda"); /* model predicted average seedling count */
	table_addcolumn("mydata","simulatedcount"); /* model output: random draw for seedling count */

	PAUSE

	return;

}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("m0",0.000010,50.0,5.0,1,0,1);
	parameter_create("beta",0.000010,50.0,5.0,1,0,1);
	parameter_create("gph",-3.0,3.0,0.10,0,0,1);
	parameter_create("gdepth",-1.0,1.0,0.10,0,0,1);
	parameter_create("gmoss",-3.0,3.0,0.10,0,0,1);

	parameter_showall();

	PAUSE

	return;
}  
/************************************************************/
void likelihood()
{
	/* writes the log-likelihood given current parameter values */
	int ii;
	double prob1;
	double lambda, ph, depth, moss, count;
	double phref=7.0, depthref=10.0, mossref=0.0;

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* loop over data */
	for(ii=0 ; ii<numdata ; ii++)
	{
		ph = table_getvalue("mydata","ph",ii);
		depth = table_getvalue("mydata","depth",ii);
		moss = table_getvalue("mydata","moss",ii);
		count = table_getvalue("mydata","count",ii);

		/* calc lambda */
		lambda = cv("beta") + cv("m0")*exp( cv("gph")*(ph-phref) + cv("gdepth")*(depth-depthref) + cv("gmoss")*(moss-mossref));
		/* prob of observed count */
		prob1 = poisson_density((int)count,lambda);

		if(prob1<0.000000010)
			prob1 = 0.000000010;

		// write prediction into data structure (nb would not be allowable under OMP parallel chains: serial chains only)
		table_writevalue("mydata","predlambda",ii,lambda);
		table_writevalue("mydata","simulatedcount",ii,poisson_draw(lambda));
		// add log prob to sum
		inc_metr_ltotnew(log(prob1));
		inc_metr_number_ok(1);
	}
	
	return;
} // end of likelihood function

/************************************************/
void final_output()
{
	/* create link to file */
	char fname[100];
		
	/* create file name for output */
	get_filzbach_path(fname, 100);
	/* now your bit */
	strcat(fname,"_my_output.txt");

	/* we're ouputting predictions */
	/* in this case, one set of preds with params set at posterior means */
	/* but could output an ensemble of predictions, with each set with param vector drawn at random */
	/* see the commented out line below... */
	params_set_to_posterior_mean();
	// params_draw_random_vector();
	/* in this case it's the likelihood function that writes the preds into the table */
	likelihood();
	/* output */
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

  fake_data(); // either fake data (to test routine) or read data (to do for real)
  read_data();
  setup_parameters();
  name_analysis("Egz2SeedrainFull");
  runmcmc(5000, 5000, 500, 500);
  final_output();
  name_analysis("Egz2SeedrainNodepth");
  parameter_fix("gdepth");
  runmcmc(5000, 50000, 5000, 5000);
  final_output();
}

#endif

