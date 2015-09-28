#include "preamble.h"
#if MODEL == 16

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
int numdata,numpp;
int link;


/************************************************************/
/* function headers                                         */
/************************************************************/

void fake_data();
void read_data();
void setup_parameters();
void fit_model();
void final_output();
double vector_probuniform(const char pname[], int vector_length, double lower_bound, double upper_bound, int sections);

/************************************************************/
void fake_data()
{
	printf("\n faking data \n");
	/* use this to create fake data to test your analysis, before you try real data */
	/* do we get these parameters back out again? */
	int pos=0,link,p1,p2,pp;
	double n1[301],n2[301]; // two dimensions this time
	double probinteracts;

	numpp=50;
	double rho1 = 0.20;
	double rho2 = 0.60;

	table_create("fakedata");
	table_addcolumn("fakedata","p1id");
	table_addcolumn("fakedata","p2id");
	table_addcolumn("fakedata","empirical");

	table_create("trueparams");
	table_addcolumn("trueparams","pid");
	table_addcolumn("trueparams","n1");
	table_addcolumn("trueparams","n2");

	// assign random parameters
	// loop over proteins
	for(pp=0 ; pp<numpp ; pp++)
	{
		table_writevalue("trueparams","pid",pp,pp);
		table_writevalue("trueparams","n1",pp,random(0.0,1.0));
		table_writevalue("trueparams","n2",pp,random(0.0,1.0));

		/* hold in array just for convenience */
		n1[pp] = table_getvalue("trueparams","n1",pp);
		n2[pp] = table_getvalue("trueparams","n2",pp);

	}

	printf("\n done random params \n");

	
	// generate fake web
	for(p1=0 ; p1<numpp-1 ; p1++)
	{
		// only do calcs for when p2>p1, avoids double counting in this symmetric web
			for(p2=p1+1 ; p2<numpp ; p2++)
			{

			// calculate prob that this ii sticks to jj
			probinteracts = 0.999 * exp( (-1.0)*pow((n1[p1]-n1[p2])/rho1,2.0)) * exp( (-1.0)*pow((n2[p1]-n2[p2])/rho2,2.0));
			// choose 1 or 0 for this link accordingly
			if(random(0.0,1.0)<probinteracts)
				link=1;
			else
				link=0;
			// write into fake data table 
			table_writevalue("fakedata","p1id",pos,p1);
			table_writevalue("fakedata","p2id",pos,p2);
			table_writevalue("fakedata","empirical",pos,link);
			// increment position (row) in table
			pos++;
		}
	}

	printf("\n done writing to table \n");

	// output fake data to file, also output true params
	table_output("trueparams","./workspace/PROTEINTrueParams.txt");
	table_output("fakedata","./workspace/PROTEINFakeData.txt");

	return;
}

/************************************************************/
void read_data()
{
	// get data from file
	table_read("mydata","./workspace/PROTEINFakeData.txt",3);
	// keep track of how many bits of data in table, i.e., how many rows
	numdata = table_numrows("mydata");
	// add column for predictions
	table_addcolumn("mydata","predprob");
	// work out number of species
	numpp = (int)table_getcolumnmax("mydata","p2id")+1; // +1 because ids run 0 ... n-1 for n species

	PAUSE

	return;
}

/************************************************************/
void setup_parameters()
{
	/* each line efines one new parameter for our model */

	/* set up n1, n2 params for proteins */
	parameter_create_vector("n1",0.0010,0.990,0.500,0,0,0,numpp);
	parameter_create_vector("n2",0.0010,0.990,0.500,0,0,0,numpp);
	/* delay n2 to aid convergence -- force to go with one dimension at first */
	parameter_delay("n2",10000);
	/* rho parameters */
	parameter_create("rho1",0.0010,10.0,1.0,1,0,1);
	parameter_create("rho2",0.0010,10.0,1.0,1,0,1);

	parameter_showall();

	PAUSE

	return;
}  
/************************************************************/
void likelihood()
{
	double probinteracts,probo1;
	int p1,p2,link;
	double ltot=0.0;
	int lcount=0;
	int ii;

	/* writes the log-likelihood given current parameter values */
	for(ii=0 ; ii<numdata ; ii++)
	{
		// get species ids for this link
		p1 = (int)table_getvalue("mydata","p1id",ii);
		p2 = (int)table_getvalue("mydata","p2id",ii);

		// calc prob of interaction according to current parameter values
		probinteracts = 0.999 * exp( (-1.0)*pow((cv("n1",p1)-cv("n1",p2))/cv("rho1"),2.0)) * exp( (-1.0)*pow((cv("n2",p1)-cv("n2",p2))/cv("rho2"),2.0));

		// write this probability into the predictions column in the mydata table
		table_writevalue("mydata","predprob",ii,probinteracts);

		// get the observed status for this link
		link = (int)table_getvalue("mydata","empirical",ii);

		// calculate probability for this observation given prediction
		if(link==1)
			probo1 = probinteracts;
		else
			probo1 = 1.0 - probinteracts;

		ltot += log(probo1);
		lcount++;

	}

	// now do bit about uniform distribution of proteins along the axes
	ltot += vector_probuniform("n1",numpp,0.0,1.0,6);
	ltot += vector_probuniform("n2",numpp,0.0,1.0,6);

	// send likelihood to MCMC algorithm
	set_metr_ltotnew(ltot);
	set_metr_number_ok(lcount);
	

	return;
}

/************************************************/
void final_output()
{	
	char fname[100];

	/* create file name for output -- outp is the 'path' */
	get_filzbach_path(fname, 100);
	/* now your bit to add to the path */
	strcat(fname,"_my_output.txt");

	/* write from internal table to external file */
	table_output("mydata",fname);

	return;

}

/****************************************************/
double vector_probuniform(const char pname[], int vector_length, double lower_bound, double upper_bound, int sections)
{
	// returns LOG probability that this vector of params has come from a uniform between lower and upper bound
	int sec, pp, sum1;
	double w = (upper_bound-lower_bound)/((double)sections);
	double expec = ((double)vector_length)/((double)sections);
	double min, max;
	double ltot=0.0;

	// loop over sections
	for(sec=0 ; sec<sections; sec++)
	{
		// set sum1 to 0
		sum1 = 0;
		// calculate bounds on section
		min = ((double)sec)*w;
		max = ((double)(sec+1))*w;
		// loop over params
		for(pp=0 ; pp<vector_length ; pp++)
		{
			if(cv(pname,pp)>=min && cv(pname,pp)<max)
				sum1++;
		}
		// do prob
		ltot += log( poisson_density(sum1,expec) );
	}

	return ltot;

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
  name_analysis("Egz6_PROTEIN");
  fake_data(); // either fake data (to test routine) or read data (to do for real)
  read_data();
  setup_parameters();
  set_chains(1);
  runmcmc(50000, 25000, 500, 500);
  final_output();
}

#endif

