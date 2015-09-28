#include "preamble.h"
#if MODEL == 13

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
	double pmax=0.60, opttemp=10.0, optprecip=500.0, sigtemp=6.0, sigprecip=500.0, occ;
	double temp, precip, probocc;
	int ii;

	/* create internal table to hold data, and model predictions (written in later) */
	table_create("fakedata");
	table_setnumrows("fakedata",299);
	table_addcolumn("fakedata","temp");
	table_addcolumn("fakedata","precip");
	table_addcolumn("fakedata","occ");
	
	/* how many data in all? */
	numdata = 299;

	for(ii=0 ; ii<numdata ; ii++)
	{
		/* draw random values for temp, precip from uniform distributions */
		temp = random(0.0,1.0)*30.0;
		precip = random(0.0,1.0)*2000.0;
		/* calculate prob occurence */
		probocc = pmax * exp((-1.0)*pow((temp-opttemp)/sigtemp,2.0));
		probocc *= exp((-1.0)*pow((precip-optprecip)/sigprecip,2.0));
		if(probocc<0.00000010)
			probocc = 0.00000010;
		/* draw 1 or 0 at random */
		if(random(0.0,1.0)<=probocc)
			occ = 1.0;
		else
			occ = 0.0;
		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata","temp",ii,temp);
		table_writevalue("fakedata","precip",ii,precip);	
		table_writevalue("fakedata","occ",ii,occ);
	}

	table_output("fakedata","./workspace/EgDataz3_SDM.txt");

	return;
}

/************************************************************/
void read_data()
{
    /* read table from external file into internal file called 'mydata' */
    table_read("mydata","./workspace/EgDataz3_SDM.txt",3);
	numdata = table_numrows("mydata");

	/* in this case, we are going to write into our internal table, the model predictions */
	/* so we need an extra column for that */
	table_addcolumn("mydata","predprob");
	table_addcolumn("mydata","predocc");

	PAUSE

	return;

}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("pmax",0.00010,0.9999,0.50,0,0,1);
	parameter_create("opttemp",0.00010,30.0,15.0,0,0,1);
    parameter_create("optprecip",0.00010,2000.0,1000.0,0,0,1);
	parameter_create("sigtemp",1.0,100.0,15.0,0,0,1);
	parameter_create("sigprecip",1.0,20000.0,1000.0,0,0,1);

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
	double probocc;
	double temp,precip,occ,randocc;

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* loop over data */
	for(ii=0 ; ii<numdata ; ii++)
	{
		/* get conditions */
		temp = table_getvalue("mydata","temp",ii);
		precip = table_getvalue("mydata","precip",ii);
		occ = table_getvalue("mydata","occ",ii);

		/* calculate prob occurence given those conditions */
		probocc = cv("pmax") * exp((-1.0)*pow((temp-cv("opttemp"))/cv("sigtemp"),2.0));
		probocc *= exp((-1.0)*pow((precip-cv("optprecip"))/cv("sigprecip"),2.0));
		if(probocc<0.00000010)
			probocc = 0.00000010;
		/* compare to real data */
		if(occ>0.0)
			prob1 = probocc;
		else
			prob1 = 1.0 - probocc;
		// add log prob to sum
		inc_metr_ltotnew(log(prob1));
		inc_metr_number_ok(1);

		// write prediction into data structure (nb would not be allowable under OMP parallel chains: serial chains only)
		table_writevalue("mydata","predprob",ii,probocc);
		// draw random occ value and write to table 
		if(random(0.0,1.0)<probocc)
			randocc=1.0;
		else
			randocc=0.0;
		table_writevalue("mydata","predocc",ii,randocc);

	}
	
	return;
}

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
  name_analysis("Egz3SDM");
  fake_data();
  read_data();
  setup_parameters();
  runmcmc(50000, 50000, 5000, 5000);
  final_output();
}

#endif

