#include "preamble.h"
#if MODEL == 4

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
	double true_lambda = 2.50;

	table_create("fakedata");
	table_addcolumn("fakedata","count");

	/* how many data in all? */
	int numdata = 99;

	for(int ii=0 ; ii<numdata ; ii++)
	{
		/* draw random value for y from normal distribution with appropriate mean and sdev */
		double count = (double)poisson_draw(true_lambda);
		table_writevalue("fakedata","count",ii,count);
	}

	// output to fake data file
	table_output("fakedata","./workspace/eg04_poisson_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read in a 1-column table from an external file, into an internal table called mydata */
	table_read("mydata","./workspace/eg04_poisson_fakedata.txt",1);

	PAUSE
	return;
}

/************************************************************/
void setup_parameters()
{	
	/* each line defines one new parameter for our model */
	/* in this case there is just one parameter in total! */
	parameter_create("lambda",0.0010,200.0,100.0,0,0,1);
	
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
		int count = (int)table_getvalue("mydata","count",ii);
		double prob1 = poisson_density(count, cv("lambda"));
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
  name_analysis("eg04_poisson_00");
  fake_data(); // either fake data (to test routine) or read data (to do for real)
  read_data();
  setup_parameters();
  //set_thinning(10);
  runmcmc(5000, 5000, 5000, 5000);
  final_output();
}

#endif

