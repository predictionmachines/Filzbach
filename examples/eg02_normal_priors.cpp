#include "preamble.h"
#if MODEL == 2

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
	/* use this function to create fake data to test your analysis, before you try real data */
	/* do we get these parameters back out again? */
	double true_sigma = 9.0; 
	double true_mean = 37.0;
	/* how many data in all? */
	int numdata = 999;

	/* create internal table to hold data */
	table_create("fakedata");
	table_addcolumn("fakedata","y");

	/* create fake data and write into internal table */
	for(int ii=0 ; ii<numdata ; ii++)
	{
		/* draw random value for y from normal distribution with appropriate mean and sdev */
		double y = normal_draw(true_mean,true_sigma);

		/* write these values to the data (as if reading real data) */
		table_writevalue("fakedata","y",ii,y);
	}

	/* output data to use later on */
	table_output("fakedata","./workspace/eg02_normal_prioirs_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read in a 1-column table from an external file, into an internal table called mydata */
    table_read("mydata","./workspace/eg02_normal_prioirs_fakedata.txt",1);

	PAUSE
	return;
}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("mean",0.0010,2000.0,100.0,0,0,1);
	parameter_create("sigma",0.0010,2000.0,20.0,1,0,1);

	/* place priors on these parameters */
	parameter_addprior("mean",500.0,100.0);
	parameter_addprior("sigma",100.0,0.10);
	
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
	for(int ii=0 ; ii<numdata ; ii++)
	{
		double y = table_getvalue("mydata","y",ii);
		double prob1 = normal_density(y, cv("mean"), cv("sigma"));
		inc_metr_ltotnew( log(prob1) );
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
  name_analysis("eg02_normal_prioirs_00");
  fake_data(); 
  read_data();
  setup_parameters();
  runmcmc(5000, 50000, 5000, 5000);
  final_output();
}

#endif
