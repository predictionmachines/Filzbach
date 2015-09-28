#include "preamble.h"
#if MODEL == 3

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
	double true_sigma1 = 9.0;
	double true_mean1 = 37.0;
	double true_sigma2 = 25.0;
	double true_mean2 = 100.0;
	double true_q = 0.80;

	/* how many data in all? */
	int numdata = 99;

	/* create internal table to hold data */
	table_create("fakedata");
	table_addcolumn("fakedata","y");

	for(int ii=0 ; ii<numdata ; ii++)
	{
		/* first, randomly choose normal1, or normal 2 */
		int choice;
		if(random(0.0,1.0)<true_q)
			choice = 1;
		else
			choice = 2;

		/* now draw random value for y from the appropriate normal distribution */
		double y;
		if(choice==1)
			y = normal_draw(true_mean1,true_sigma1);
		if(choice==2)
			y = normal_draw(true_mean2,true_sigma2);

		/* write the y value to the data (as if reading real data) */
		table_writevalue("fakedata","y",ii,y);
	}

	/* output data to use later on */
	table_output("fakedata","./workspace/eg03_mixednormal_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read from external file to internal table */
    table_read("mydata","./workspace/eg03_mixednormal_fakedata.txt",1);

	PAUSE
	return;
}

/************************************************************/
void setup_parameters()
{
	/* each line efines one new parameter for our model */

	parameter_create("mean1",0.0010,2000.0,100.0,0,0,1);
	parameter_create("sigma1",0.0010,2000.0,20.0,1,0,1);
	parameter_create("mean2",0.0010,2000.0,100.0,0,0,1);
	parameter_create("sigma2",0.0010,2000.0,20.0,1,0,1);
	parameter_create("q",0.00010,0.999,0.50,0,0,1);

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
		/* if this val came from dist 1 */
		double prob1 = normal_density(y, cv("mean1"), cv("sigma1")); 
		/* if this val came from dist 1 */
		double prob2 = normal_density(y, cv("mean2"), cv("sigma2"));
		/* weighted prob */
		double prob3 = prob1*cv("q") + prob2*( 1.0 - cv("q") );
		inc_metr_ltotnew(log(prob3));
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
  name_analysis("eg03_mixednormal_00");
  fake_data();
  read_data();
  setup_parameters();
  runmcmc(5000, 50000, 5000, 5000);
  final_output();
}

#endif

