#include "preamble.h"
#if MODEL == 0

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "filzbach.h"

//#define PAUSE // uncomment to run without pauses
void pause(){PAUSE}


/************************************************************/
void fake_data()
{
	/* use this function to create fake data to test your analysis, before you try real data */

	/* how many data in all? */
	int numdata = 999;

	/* create internal table to hold data */
	table_create("fakedata");

	/* define table structure */ 
	// table_addcolumn("fakedata", <column_name>);
	// ...

	/* create fake data */
	for(int ii = 0; ii < numdata; ii++)
	{
		/* compute values of a fake data row */

		/* write these values to the data (as if reading real data) */
		// table_writevalue("fakedata", <column_name>, ii, <computed_value>);
		// ...
	}

	/* output data to use later on */
	table_output("fakedata","./workspace/mymodel_fakedata.txt");

	return;
}

/************************************************************/
void read_data()
{
	/* read in a 1-column table from an external file, into an internal table called mydata */
	/* change the second parameter when reading real data instead of fake data */
	// TODO: change the third parameter to the actual number of columns!
	table_read("mydata", "./workspace/mymodel_fakedata.txt", 0);

	/* add columns for the likelihood function to write predictions */ 
	// table_addcolumn("mydata", <column_name>);
	// ...

	PAUSE // examine printout to visually check the proper data was read
	return;
}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	// parameter_create(<param_name>,<lower_bound>,<upper_bound>,<initial_value>,<type>,<is_fixed?>,<is displayed>);

	parameter_showall();

	PAUSE // examine printout to visually check parameters as seen by Filzbach
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
		/* read data stored in a table */
		// double <data_var> = table_getvalue("mydata",<column_name>,ii);

		/* read current parameter values */
		// double <param_var> = cv(<param_name>);

		/* you may compute and write some "predictions" to compare the data with */
		// double <pred_var> = <expression_of_data_and_param>;
		// table_writevalue("mydata", <coulumn_name>, ii, <pred_var>

		/* now we have all we need to compute the likelihood of this data row */
		// double prob = <expression_of_data_param_and_pred>;

		/* increment overall log likelihood and sample size */
		//inc_metr_ltotnew( log(prob) );
		//inc_metr_number_ok(1);
	}

	return;
}

/************************************************/
void final_output()
{
	char fname[100];

	/* create file name for output -- fname gets the 'path' and analysis name */
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
	atexit(pause); // sets useful pause when debugging in Visual Studio
	// set the likelihood function pointer
	pfn_likelihood = &likelihood;

	initialize_filzbach();
	/* analysis name becomes part of all file names written by Filzbach */
	name_analysis("mymodel_00");
	fake_data(); 
	read_data();
	setup_parameters();
	set_chains(3);
	// runmcmc(10000, 10000, 10000, 10000);
	final_output();
}

#endif
