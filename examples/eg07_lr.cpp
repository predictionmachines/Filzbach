#include "preamble.h"
#if MODEL == 7

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
	double true_c = 25.0;
	double true_m = 15.0;
	double true_sigma = 100.0;

	/* how many data? */
	int numdata = 99;

	/* create internal table to hold data. and model predictions (written in later) */
	table_create("fakedata");
	table_addcolumn("fakedata","x");
	table_addcolumn("fakedata","y");

	for(int ii = 0; ii < numdata; ii++)
	{
		/* draw random value for x from uniform distribution */
		double x = random(0.0,100.0);
		/* calculte the mean y for this value of x */
		double mean = true_c + true_m*x;
		/* draw a random y from a normal with this mean, and with sdev sigma */
		double y = normal_draw(mean, true_sigma);
		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata", "x", ii, x);
		table_writevalue("fakedata", "y", ii, y);		
	}

	/* output data to use later on */
	table_output("fakedata","./workspace/eg07_lr_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read table from external file into internal file called 'mydata' */
	table_read("mydata","./workspace/eg07_lr_fakedata.txt",2);

	/* in this case, we are going to write into our internal table, the model predictions */
	/* so we need an extra columnd for that */
	table_addcolumn("mydata","predy1");
	table_addcolumn("mydata","predy2");

	PAUSE
		return;
}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("c",-2000.0,2000.0,1000.0,0,0,1);
	parameter_create("m",1.0,200.0,10.0,0,0,1);
	parameter_create("sigma",0.0010,200.0,20.0,1,0,1);

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
		double predy = cv("c")+cv("m")*x;
		double prob1 = normal_density(y, predy, cv("sigma"));
		// write prediction into data structure (nb would not be allowable under OMP parallel chains: serial chains only)
		table_writevalue("mydata", "predy1", ii, predy); // expectation
		table_writevalue("mydata", "predy2", ii, normal_draw(predy,cv("sigma"))); // random draw including unexplained variation
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
	name_analysis("eg07_lr_00");
	fake_data();
	read_data();
	setup_parameters();
	set_chains(3);
	runmcmc(5000, 50000, 5000, 5000);
	final_output();
}

#endif

