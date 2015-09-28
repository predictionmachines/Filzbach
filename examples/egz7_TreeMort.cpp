#include "preamble.h"
#if MODEL == 17

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
	double k0 = -4.0, k1 = 0.20, pdead, mu, dead, tinterval, dbh;
	int ii;

	/* create internal table to hold data. and model predictions (written in later) */
	table_create("fakedata");
	table_addcolumn("fakedata","dbh");
	table_addcolumn("fakedata","tinterval");
	table_addcolumn("fakedata","dead");

	/* how many data? */
	int numdata = 500;

	for(ii=0 ; ii<numdata ; ii++)
	{
		dbh = random(0.0,50.0);
		tinterval = random(1.0,12.0);
		/* calculate annual prob of death */
		mu = logistic(k0 + k1*dbh);
		// overall prob of death
		pdead = 1.0 - pow(1.0 - mu,tinterval);
		if(random(0.0,1.0)<pdead)
			dead=1.0;
		else
			dead=0.0;
		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata","dbh",ii,dbh);
		table_writevalue("fakedata","tinterval",ii,tinterval);	
		table_writevalue("fakedata","dead",ii,dead);	
	}

    /* read table from external file into internal file called 'mydata' */
	char fname[100];
	get_filzbach_path(fname, 100);
	/* now your bit */
	strcat(fname,"_my_fake_data.txt");
	table_output("fakedata",fname);

	return;
}

/************************************************************/
void read_data()
{
    /* read table from external file into internal file called 'mydata' */
	char fname[100];
	get_filzbach_path(fname, 100);
	/* now your bit */
	strcat(fname,"_my_fake_data.txt");
    table_read("mydata", fname, 3); // number of columns -- sorry!
	int numdata = table_numrows("mydata");

	/* in this case, we are going to write into our internal table, the model predictions */
	/* so we need an extra columnd for that */
	table_addcolumn("mydata","pdead");

	PAUSE

	return;

}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("k0",		-5.0, 5.0,	-3.0, 0, 0, 1);
	parameter_create("k1",		-0.50, 0.50, 0.0, 0,	0,	1);

	parameter_showall();

	PAUSE

	return;
}  
/************************************************************/
void likelihood()
{
	/* writes the log-likelihood given current parameter values */
	double mu, prob1;
	double dbh, tinterval, pdead, dead;

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* loop over data */
	int numdata = table_numrows("mydata");
	for(int ii=0 ; ii<numdata ; ii++) // equiv. to sigma sign i=1 to n
	{
		dbh = table_getvalue("mydata","dbh",ii);
		tinterval = table_getvalue("mydata","tinterval",ii);
		dead = table_getvalue("mydata","dead",ii);

		mu = logistic( cv("k0") + cv("k1")*dbh);
		pdead = 1.0 - pow(1.0-mu,tinterval);
		if(dead>0.0)
			prob1 = pdead;
		else
			prob1 = 1- pdead;

		// write prediction into data structure (nb would not be allowable under OMP parallel chains: serial chains only)
		table_writevalue("mydata","pdead",ii,pdead); // expectation

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

	/* create file name for output */
	get_filzbach_path(fname, 100);
	/* now your bit */
	strcat(fname,"_my_output.txt");
	
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
	name_analysis("egz7_TreeMort");
	fake_data();
	read_data();
	setup_parameters();
	runmcmc(50000, 50000, 50000, 50000);
	final_output();
}

#endif

