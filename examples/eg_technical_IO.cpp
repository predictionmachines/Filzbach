#include "examples.h"
#ifdef TECHNICAL_IO

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 #include <time.h>;
#include <string.h>

#include "filzbachtbl.h";
#include "filzbachstat.h";
#include "filzbachmetr.h";

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
void setup_model();
void fit_model();
void final_output();

/************************************************************/
void fake_data()
{	
	table *ptrtable;
	table mytable;

	int tid, cc, ii;

	double temp1,temp2;
	double adouble = 7.99;
	int anint = 2;
	int tspp;

	/* create two different internal tables */
	table_create("countdata");
	table_setnumrows("countdata",512); // name, and number of rows
	table_addcolumn("countdata","species");
	table_addcolumn("countdata","count");

	table_create("heightdbhdata");
	table_setnumrows("heightdbhdata",2160); // name and number of rows are different
	table_addcolumn("heightdbhdata","species");
	table_addcolumn("heightdbhdata","dbh");
	table_addcolumn("heightdbhdata","height");

	/* write a value into the first table */
	table_writevalue("countdata","species",15,126); // writes '126' into the 15th row, of the column called 'species' in the table called 'coutdata' */
	temp1 = 126;
	/* read this value back in */
	temp2 = table_getvalue("countdata","species",15);
	printf("\n testing using tables: these two numbers should be the same: %lf, %lf",temp1,temp2);
	
	/* output the second table to a file */
	table_output("heightdbhdata","./workspace/egtechnical_testoutput1.txt");

	/* read a THIRD table from the file we just made... */
	table_read("heightdbhdata_copy","./workspace/egtechnical_testoutput1.txt",3); // 3 is the no. columns

	/* output this to a different file! */
	table_output("heightdbhdata_copy","./workspace/egtechnical_testoutput1_copy.txt");

	/* write '1' into the species column of every row in countdata -- the slow, easy way */
	for(ii=1 ; ii<=512 ; ii++)
	{
		table_writevalue("countdata","species",ii,1.0);
	}
		
	/* write '1' into the species column of every row in countdata -- the harder, faster way */
	tid = table_getID("countdata");
	cc = table_getcolumnID(tid,"species");
	ptrtable = table_get(tid);
	mytable = *ptrtable;

	for(ii=1 ; ii<=512 ; ii++)
	{		
		mytable.data[cc][ii] = 1.0;
	}

	/* read contents of the same column -- the hard, fast way */
	tid = table_getID("countdata");
	cc = table_getcolumnID(tid,"species");
	ptrtable = table_get(tid);
	mytable = *ptrtable;

	for(ii=1 ; ii<=512 ; ii++)
	{			
		tspp = mytable.data[cc][ii];
	}
	/* print some numbers to screen */
	printf("\n here is a double: %lf, here is an integer: %d \n",adouble, anint);

	return;
}

/************************************************************/
void read_data()
{
	table_read("mydata","./workspace/EgData4_poisson.txt",1);
	numdata = table_numrows("mydata");
	printf("\n read data OK, read %d samples \n \n",numdata);
	system ("pause");

	return;
}

/************************************************************/
void setup_model()
{
	int mm;
	/* each line defines one new parameter for our model */
	parameter_create_vector("mean",0.0010,200.0,100.0,1,0,1,10);
	
	parameter_showall();

	system("pause");

	return;
}  
/************************************************************/
void lnlike()
{
	/* writes the log-likelihood given current parameter values */
	int ii,spp;
	double prob1,y;
	double mean, sigma;

	/* set sum over log-likelihood to zero */
	set_metr_ltotnew(0.0);
	set_metr_number_ok(0);

	/* get model parameters from list held by metropolis header */
	

	/* loop over data */
	for(ii=1 ; ii<=numdata ; ii++)
	{
		/* get observed count and species id */
		y = table_getvalue("mydata","y",ii);
		spp = (int)table_getvalue("mydata","species",ii);
		/* get parameter for this species */
		mean = parameter_getvalue("mean",spp);
		prob1 = poisson_density(y, mean);
		inc_metr_ltotnew(log(prob1));
		inc_metr_number_ok(1);
	}
	
	return;
} // end of lnlike function

/************************************************/
void final_output()
{
	char fname[100];
		
	/* create file name for output */
	/* automatic bit */
	strcpy(fname,outp);
	/* now your bit */
	strcat(fname,"_my_output.txt");

	/* output */
	table_output("mydata",fname);

	return;

}

/**************************************************************/
int main()
/* the control function that calls everything else */
{
  // initialize_filzbach();
  // name_analysis("Eg4_2_PoissonMultispp");
  fake_data(); // either fake data (to test routine) or read data (to do for real)
  // read_data();
  // setup_model();
  // runmcmc(5000, 5000, 500, 500);
  // final_output();

  system ("pause");
}

#endif
