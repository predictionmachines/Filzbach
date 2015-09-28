#include "preamble.h"
#if MODEL == 11

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

/* custom functions */
double model_biomass(int ii, int code);

/************************************************************/
void fake_data()
{
	/* use this to create fake data to test your analysis, before you try real data */
	/* do we get these parameters back out again? */
	double meanmassabv, meanmassblw, biomassabv, biomassblw, G0=0.05, vinside=0.20, vfert=0.40, vshade=-0.30, mref=5.0, sigma = 0.20;
	double potvolume,alpha=5.0,beta=0.50;
	double inside, fert, shade, day;
	double sigmaabv,sigmablw;
	int ii;

	/* create internal table to hold data */
	table_create("mydata");
	table_setnumrows("mydata", 500);
	table_addcolumn("mydata", "inside");
	table_addcolumn("mydata", "fert");
	table_addcolumn("mydata", "shade");
	table_addcolumn("mydata", "potvolume");
	table_addcolumn("mydata", "day");
	table_addcolumn("mydata", "biomassabv");
	table_addcolumn("mydata", "biomassblw");
	table_addcolumn("mydata", "predmass1"); // this for model predictions
	table_addcolumn("mydata", "predmass2"); // this for model predictions

	sigmaabv = parameter_getvalue("sigmaabv");
	sigmablw = parameter_getvalue("sigmablw");
	
	for(ii=1 ; ii<=500 ; ii++)
	{
	    /* draw random conditions and sample day */
		if(random(0.0,1.0)<0.50)
			inside = 1.0;
		else
			inside = 0.0;
		if(random(0.0,1.0)<0.50)
			shade = 1.0;
		else
			shade = 0.0;
		potvolume = ((double)((int)(random(0.0,1.0)*5.0)))/2.0;
		fert = (double)((int)(random(0.0,1.0)*5.0));
		day = random(0.0,1.0)*100.0;

		/* write the y value to the data (as if reading real data) */
		table_writevalue("mydata","inside",ii,inside);
		table_writevalue("mydata","shade",ii,shade);
		table_writevalue("mydata","fert",ii,fert);
		table_writevalue("mydata","potvolume",ii,potvolume);
		table_writevalue("mydata","day",ii,day);

		/* get expected size after correct no. days */
		meanmassabv = model_biomass(ii,0);
		meanmassblw = model_biomass(ii,1);

		/* generate actual biomass as a draw from a log-normal around this expectation */
		biomassabv = exp(normal_draw(log(meanmassabv),sigmaabv));
		biomassblw = exp(normal_draw(log(meanmassblw),sigmablw));

		/* write fake biomasses into data table */
		table_writevalue("mydata","biomassabv",ii,biomassabv);
		table_writevalue("mydata","biomassblw",ii,biomassblw);
	}

	/* how many data in all? */
	numdata = 500;

	/* immediately output fake data */
	table_output("mydata","./workspace/fake_data.txt");

	return;
}

/************************************************************/
double model_biomass(int ii, int code)
{
	int dd=0;
	double grate, grateabv, grateblw, massabv, massblw;
	double f;
	double imassabv, imassblw, G0, vinside, vshade, vfert, alpha, beta, gamma0;
	double inside, shade, fert, potv, day, gmax, mref;
	double gamma1;

	/* get model parameters from list held by metropolis header */
	imassabv = parameter_getvalue("inmassabv");
	imassblw = parameter_getvalue("inmassblw");
	G0 = parameter_getvalue("G0");
	// mref = parameter_getvalue("mref");
	vinside = parameter_getvalue("vinside");
	vshade = parameter_getvalue("vshade");
	vfert = parameter_getvalue("vfert");
	alpha = parameter_getvalue("alpha");
	beta = parameter_getvalue("beta");
	gamma0 = parameter_getvalue("gamma0");
	gamma1 = parameter_getvalue("gamma1");

	/* get conditions for this plant */
	inside = table_getvalue("mydata","inside",ii);
	shade = table_getvalue("mydata","shade",ii);
	fert = table_getvalue("mydata","fert",ii);
	day = table_getvalue("mydata","day",ii);
	potv = table_getvalue("mydata","potvolume",ii);

	/* calculate gmax */
	gmax = G0 * exp( vinside*inside + vshade*shade + vfert*fert );
	gmax *= pow(potv,beta);
	/* calculate mrepow(potv,beta);f as a function of pot volume */
	mref = alpha * pow(potv,beta);

	massabv = imassabv;
	massblw = imassblw;
	for(dd=0 ; dd<(int)day ; dd++)
	{
		if(massabv<mref)
			grate = gmax * (massabv/mref);
		else
			grate = gmax;

		f = logistic(gamma0 + gamma1*fert);

		grateabv = f * grate;
		grateblw = (1.0-f) * grate;

		massabv+=grateabv;
		massblw+=grateblw;
	}

	if(code==0)
		return massabv;
	else
		return massblw;
}

/************************************************************/
void read_data()
{
	/* read from external file to internal table */

    table_read("mydata","./workspace/EgDataz1_plantgrowth.txt",5);
	numdata = table_numrows("mydata");

	table_addcolumn("mydata","predmass1"); // this for model predictions
	table_addcolumn("mydata","predmass2"); // this for model predictions


	return;

}

/************************************************************/
void setup_parameters()
{
	/* each line defines one new parameter for our model */

	parameter_create("inmassabv",  0.00010, 1.0,  0.25,  1, 0, 1);
	parameter_create("inmassblw",  0.00010, 1.0,  0.25,  1, 0, 1);
	parameter_create("G0",         0.00010, 1.0,  0.050, 1, 0, 1);
//  parameter_create("mref",       0.10,    20.0, 2.0,   1, 0, 1);
	parameter_create("vinside",   -2.0,     2.0,  0.20,  0, 0, 1);
	parameter_create("vshade",    -2.0,     2.0, -0.20,  0, 0, 1);
	parameter_create("vfert",     -2.0,     2.0,  0.30,  0, 0, 1);
	parameter_create("sigmaabv",   0.0010,  2.0,  0.20,  1, 0, 1);
	parameter_create("sigmablw",   0.0010,  2.0,  0.20,  1, 0, 1);
	parameter_create("alpha",      0.10,    20.0, 2.0,   1, 0, 1);
	parameter_create("beta",       0.10,    1.90, 0.50,  1, 0, 1);
	parameter_create("gamma0",    -5.0,     5.0,  1.0,   0, 0, 1);
	parameter_create("gamma1",    -0.50,    0.50, 0.10,  0, 0, 1);

	parameter_showall();

	PAUSE

	return;
}  
/************************************************************/

void likelihood()
{
	/* writes the log-likelihood given current parameter values */
	int ii;
	double prob1,prob2;
	double omassabv, omassblw;
	double sigmaabv, sigmablw;
	double pmassabv, pmassblw;

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* get current values for noise parameters */
	sigmaabv = parameter_getvalue("sigmaabv");
	sigmablw = parameter_getvalue("sigmablw");

	/* loop over data */
	for(ii=1 ; ii<=numdata ; ii++)
	{

		/* get observed biomass */
		omassabv = table_getvalue("mydata","biomassabv",ii);
		omassblw = table_getvalue("mydata","biomassblw",ii);


		/* predict expected final biomass */
		pmassabv = model_biomass(ii, 0);
		pmassblw = model_biomass(ii, 1);

		/* write pmass into prediction column (nb would not be allowable under OMP parallel chains: serial chains only) */
		table_writevalue("mydata","predmass1",ii,pmassabv);
		/* write the second column, which has noise like the real data */
		table_writevalue("mydata","predmass2",ii,pmassblw);

		/* get probability for these observed masses */
		prob1 = normal_density(log(omassabv),log(pmassabv),sigmaabv);
		prob2 = normal_density(log(omassblw),log(pmassblw),sigmablw);
		inc_metr_ltotnew(log(prob1) + log(prob2));
		inc_metr_number_ok(2);
	}
	
	return;
}

/************************************************/
void final_output()
{
	char fname[100];
		
        /* create file name for output -- outp is the 'path' */
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
  name_analysis("Egz1PlantGrowth");
  setup_parameters();
  fake_data(); // either fake data (to test routine) or read data (to do for real)
  // read_data();
  runmcmc(50000, 50000, 5000, 5000);
  final_output();
}

#endif

