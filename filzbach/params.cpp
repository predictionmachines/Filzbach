#include "stdafx.h"
#include "filzbach.h"
#include "dwp_metropolis.h"
#include "params.h"

// the reference paramter chain
param *base_para = NULL;
// the parameter chains
param **mpara = NULL;
// number of paramters
int paramcount = 0;

/* 
 * Returns the current paramter set depending whether chains are computed or not.
 * If chains are computed, the corresponding paramter set is returned 
 * otherwise the paramters of the base chain.
 */
param* current_para(int currentchain)
{
	if (currentchain < 0)
	{
		return base_para;
	}
	else
	{
		return mpara[currentchain];
	}	
}

param* current_para()
{
	return current_para(get_currentchain());
}

/*
* Allocates memory for an additional parameter.
*/
void palloc()
{
	if (paramcount == 0)
	{
		// initialize base chains
		delete[] base_para;
		base_para = new param[1];
		CHECK(base_para!=NULL, "Memory allocation failed.");
	}
	else 
	{
		// allocate an additioanl paramter for the base chain
		param *tempbase_para = new param[paramcount + 1];
		CHECK(tempbase_para!=NULL, "Memory allocation failed.");
		for (int p = 0; p < paramcount; p++)
		{
			tempbase_para[p] = base_para[p];
		}

		delete [] base_para;
		base_para = tempbase_para;	
	}

	paramcount++;
}

char printbuf[256];


/// <summary>Display name of a parameter</summary>
/// <returns>An allocated C-string</returns>
/// <remarks>Call delete[] ret-value; when it is not needed any more.</remarks>
char* pprintname(int index)
{
	param* params = base_para;
	ASSERT(params!=NULL, "params array not specified");
	ASSERT((index>=0) && (index<paramcount), "index out of range");
	int i = 1;
	while ((index-i>=0) && (strcmp(params[index].name, params[index-i].name)==0)) i++;
	if (i>1 || ((index+1<paramcount) && (strcmp(params[index].name, params[index+1].name)==0)))
	{
		sprintf_s(printbuf, "%s[%d]", params[index].name, i-1);
		size_t len = strlen(printbuf)+1;
		char* result = new char[len];
		strcpy_s(result, len, printbuf);
		return result;
	}
	else
	{
		size_t len = strlen(params[index].name)+1;
		char* result = new char[len];
		strcpy_s(result, len, params[index].name);
		return result;
	}

}

/*
* Initializes a paramter with the given values.
* The initialization affects onlz the base chain.
*/
void metr_setup_param(param* p, const char name[], double lb, double ub, double value, int type, int fixed, int delay, int display)
{	
	strcpy(p->name, name);

	p->lb = lb;
	p->ub = ub;
	p->value = value;
	p->fixed = fixed;
	p->type = type;
	p->delay = delay;

	if(p->delay == 1)
	{
		p->delay = DELAY;
	}

	p->display = display;

	p->priormean = 0.0;
	p->priorsdev=0.0;
	p->prioryn=0;
	p->rootRhat=-999.0;
}

/*
* Creates a new paramter with the given values.
* The given paramter is created for all chains.
*/
int parameter_create(const char name[], double lb, double ub, double val, int type, int fixed, int dsply)
{
	CHECK(name!=NULL, "Parameter name cannot be null.");
	CHECK(ub>lb, "Upper bound must be greater than lower bound.");
	CHECK(type<=0 || (ub>0 && lb > 0), "For type 1 parameters both lower and upper bounds must be positive.");
	palloc();
	
	metr_setup_param(&base_para[paramcount - 1], name, lb, ub, val, type, fixed, 0, dsply);

	base_para[paramcount - 1].ival = val;

	return paramcount - 1;
}


/*
* Create an vector of parameter called 'name', all with same info (bounds etc.). 
*/
void parameter_create_vector(const char name[], double lb, double ub, double val, int type, int fixed, int dsply, int number)
{		
	for(int n = 0 ; n < number ; n++)
	{
		parameter_create(name, lb, ub, val, type, fixed, dsply);
	}
}

/***************************************************************************/
double cv(const char name[])
{
	int hits=0;
	int hitmm=0;

	for(int mm=0 ; mm < paramcount ; mm++)
	{
		if(strcmp(current_para()[mm].name,name)==0)
		{
			hits++;
			hitmm = mm;
		}
	}

	if (hits==1)
	{
		return current_para()[hitmm].value;
	}
	else if (hits==0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
	}
	else if (hits>1)
	{
		sprintf_s(printbuf, "More than than one parameter named %s. Did you mean a vector of parameter?", name);
	}
	CHECK(false, printbuf);

}

/************************************************************************/
double cv(int mm)
{
	if(mm>=0 && mm<paramcount)
	{
		return current_para()[mm].value;
	}
	else
	{
		sprintf_s(printbuf, "error with CV, no param numbered %d \n",mm);
		CHECK(false, printbuf);
	}	
}

/************************************************************************/
double parameter_getvalue(const char name[])
{
	return cv(name);
}


/********************************************************/
void parameter_fix(const char name[])
{
	int mm,hits=0;
	
	param *params = current_para();
	for(mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			params[mm].fixed=1;
			params[mm].value = params[mm].ival;
		}
	}

	if(hits<1)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
}

/********************************************************/
void parameter_fix(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		current_para()[mm].fixed = 1;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/********************************************************/
void parameter_fix_value(const char name[], double val)
{
	int mm,hits=0;

	param *params = current_para();
	for(mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			params[mm].fixed=1;
			params[mm].value = val;
		}
	}

	if(hits<1)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
}

/********************************************************/
void parameter_fix_value(const char name[], int number, double val)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		current_para()[mm].fixed = 1;
		current_para()[mm].value = val;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/********************************************************/
void parameter_setvalue(const char name[], double val)
{
	param *params = current_para();
	for(int mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			params[mm].value = val;
			return;
		}
	}
	sprintf_s(printbuf, "No parameter named %s found.", name);
	CHECK(false, printbuf);
}

/********************************************************/
void parameter_setvalue(int mm, double val)
{
	if(mm<=0 && mm<paramcount)
	{
		param *params = current_para();
		params[mm].value = val;
		return;
	}
	else
	{
		sprintf_s(printbuf, "error with setvalue, no param numbered %d",mm);
		CHECK(false, printbuf);
	}

}

/********************************************************/
void parameter_setvalue(const char name[], int index, double val)
{
	int mm,hits=index;

	param *params = current_para();
	for(mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
			if (hits-- == 0)
			{
				params[mm].value = val;
				return;
			}
	}
	sprintf_s(printbuf, "error, parameter setvalue, no parameter named %s with index %d found.", name, index);
	CHECK(false, printbuf);
}

/********************************************************/
void parameter_unfix(const char name[])
{
	int mm,hits=0;

	param *params = current_para();
	for(mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			params[mm].fixed=0;
		}
	}

	if(hits<1)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
}

/********************************************************/
void parameter_unfix(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		current_para()[mm].fixed = 0;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
void parameter_delay(const char name[], int iters)
{
	int hits=0;

	param *params = current_para();
	for(int mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;			
			params[mm].delay=iters;
		}
	}

	if(hits<1)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
int parameter_ifchanged(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else{
		if(current_para()[ii].alt>0 || current_para()[ii].decis==1){ // if changed this time, or changed and reject last time
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

/************************************************************************/
double l95(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bcred_l95;
	}
}

/************************************************************************/
double l95(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bcred_l95;
}

/************************************************************************/
double l68(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bcred_l68;
	}
}

/************************************************************************/
double l68(int ii)
{

	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bcred_l68;
}

/************************************************************************/
double postmedian(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bcred_median;
	}
}

/************************************************************************/
double postmedian(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bcred_median;
}

/************************************************************************/
double u68(const char name[])
{		
	int ii = parameter_returnindex(name);

	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bcred_u68;
	}
}

/************************************************************************/
double u68(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bcred_u68;
}

/************************************************************************/
double u95(const char name[])
{
	int ii = parameter_returnindex(name);

	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bcred_u95;
	}
}

/********* Profile likelihood ************/
double pl_l95(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].MLEl95;
	}
}

double pl_l68(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].MLEl68;
	}
}

double mle(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].MLE;
	}
}

double pl_u68(const char name[])
{		
	int ii = parameter_returnindex(name);

	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].MLEu68;
	}
}

double pl_u95(const char name[])
{
	int ii = parameter_returnindex(name);

	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].MLEu95;
	}
}


/************************************************************************/
double u95(int ii)
{

	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bcred_u95;
}

/************************************************************************/
double postmean(const char name[])
{
	int ii;

	ii = parameter_returnindex(name);
	if(ii<0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	else
	{
		return current_para()[ii].bayes_mean;
	}
}

/************************************************************************/
double postmean(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].bayes_mean;
}

/************************************************************************/
double l95(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return l95(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double l68(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return l68(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double postmedian(const char name[], int number)
{	
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return postmedian(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double u68(const char name[], int number)
{	
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return u68(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double u95(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return u95(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double postmean(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return postmean(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}


/********************* Profile likelihood results ***************************/
double pl_l95(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].MLEl95;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

double pl_l68(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].MLEl68;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

double mle(const char name[], int number)
{	
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].MLE;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

double pl_u68(const char name[], int number)
{	
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].MLEu68;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

double pl_u95(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].MLEu95;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
int parameter_ifchanged(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	if(current_para()[ii].alt > 0 || current_para()[ii].decis==1)
	{ 
		return 1;
	}
	else
	{
		return 0;
	}
}

/************************************************************************/
int parameter_ifchanged(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return parameter_ifchanged(mm);
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double cv(const char name[], int number)
{
	int mm = parameter_returnindex(name, number);

	if(mm>=0)
	{
		return current_para()[mm].value;
	}
	else
	{
		sprintf_s(printbuf, "No parameter named %s found or it is shorter than %d.", name, number);
		CHECK(false, printbuf);
	}
}

/************************************************************************/
double parameter_getvalue(const char name[], int number)
{
	return cv(name,number);
}

/*
 * Find the 'number'th copy of a parameter called 'name' and report the index.
 */
int parameter_returnindex(const char name[], int number)
{	
	param* params = current_para();
	for (int p = 0; p < paramcount; ++p)
	{
		if(strcmp(params[p].name, name) == 0)
		{
			if (number==0) return p;
			else if (strcmp(params[p+number].name, name) == 0) return p+number;
			else return -1;
		}
	}
	return -1;
}

/***************************************************************************/
double parameter_getvalue(int ii)
{
	CHECK(ii>=0 && ii<paramcount,"Index out of range");
	return current_para()[ii].value;
}

/***************************************************************************/
int parameter_returnindex(const char name[])
{
	int hits=0,hitmm=-1;

	param* params = current_para();
	for(int mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			hitmm = mm;
		}
	}

	if(hits==1)
	{
		return hitmm;
	}

	return -1;
}

/*
 * Outputs parameters (including details like lb, ub) to screen.
 * Must be called from the likelihood function only.
 */
void parameter_showall()
{
	printf("\nshowing all parameters: \n");

	param* p = current_para();
	for(int mm=0 ; mm<paramcount ; mm++)
	{
		char* print_name = pprintname(mm);
		printf("%s \t lb,value,ub:(%lf, %lf, %lf)\n", print_name, p[mm].lb, p[mm].value, p[mm].ub);
		printf("\t type,fixed,dsply:(%d, %d, %d)\n", p[mm].type, p[mm].fixed, p[mm].display);
		printf("\t prior:1=yes(%d) mean,sdev(%lf,%lf)\n", p[mm].prioryn, p[mm].priormean, p[mm].priorsdev);
		delete print_name;
	}	
}

/***************************************************************************/
void parameter_addprior(const char name[], double mean, double sdev)
{	
	int hits = 0;
	int hitmm = 0;

	for(int mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(current_para()[mm].name,name)==0)
		{
			hits++;
			hitmm = mm;
		}
	}

	if(hits==1)
	{
		current_para()[hitmm].prioryn = 1;
		current_para()[hitmm].priormean = mean;
		current_para()[hitmm].priorsdev = sdev;
		return;
	}
	else if (hits==0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
	}
	else if (hits>1)
	{
		sprintf_s(printbuf, "More than than one parameter named %s. Did you mean a vector of parameter?", name);
	}
	CHECK(false, printbuf);
}

/***************************************************************************/
void parameter_addprior_vector(const char name[], double mean, double sdev)
{
	int mm,hits=0,hitmm=0;

	param* params = current_para();
	for(mm = 0; mm < paramcount; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			params[hitmm].prioryn = 1;
			params[hitmm].priormean = mean;
			params[hitmm].priorsdev = sdev;			
		}
	}

	if(hits==0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
}

/***************************************************************************/
void parameter_normalize_vector(const char name[], double mean, double sdev)
{
	int mm,hits;
	double sum1=0.0,sum2=0.0,pdev=0.0;

	hits=0;
	param* params = current_para();
	for(mm=0; mm<paramcount; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			hits++;
			sum1 += params[mm].value;
			sum2 += 1.0;
		}
	}

	if(hits==0)
	{
		sprintf_s(printbuf, "No parameter named %s found.", name);
		CHECK(false, printbuf);
	}
	if(hits==1)
	{
		sprintf_s(printbuf, "The parameter named %s is a scalar.", name);
		CHECK(false, printbuf);
	}

	/* correct mean */
	sum1 /= sum2;
	for(mm=0; mm < paramcount; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			params[mm].value += mean - sum1;
			pdev += pow(params[mm].value-mean,2.0);
		}
	}
	pdev/=sum2-1.0;
	pdev=sqrt(pdev);

	// if pdev ridiculously small bail now
	if(pdev<0.00010*mean)
		return;

	/* correct pdev */
	for(mm=0 ; mm < paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			params[mm].value = mean + (params[mm].value-mean)*(sdev/pdev);
		}
	}

	/* check for straying outside bounds */
	for(mm=0 ; mm<paramcount ; mm++)
	{
		if(strcmp(params[mm].name,name)==0)
		{
			if(params[mm].value<params[mm].lb)
				params[mm].value = params[mm].lb;
			if(params[mm].value>params[mm].ub)
				params[mm].value = params[mm].ub;
		}
	}

	return;
}

/***************************************************************************/
void parameter_addprior(int mm, double mean, double sdev)
{
	param* p = current_para()+mm;
	CHECK(mm>=0 && mm<paramcount,"Index out of range");
	p->prioryn = 1;
	p->priormean = mean;
	p->priorsdev = sdev;
	return;
}

/* 
* Allocates and initializes paramter vectors for the number or required chain.
* If nothing was specified by the user, one single chain is instantiated.
*/
void init_chains(int chaincount)
{	
	ASSERT(chaincount>0,"chain count not positive");

	static int mpara_chaincount = 0;
	if (mpara != NULL) for (int c=0; c<mpara_chaincount; c++) delete[] mpara[c];
	delete[] mpara;
	mpara = new param *[chaincount];
	mpara_chaincount = chaincount;
	CHECK(mpara!=NULL,"Memory allocation failed.");

	// for all chains
	for (int c = 0; c < chaincount; c++)
	{
		// create new chain
		mpara[c] = new param[paramcount];
		CHECK(mpara[c]!=NULL,"Memory allocation failed.");

		// initialize the paramters
		//pinitialize(c, paramcount - 1); //VL it doesn't do the job anyway

		// copy all prameters from the base chain
		for (int p = 0; p < paramcount; p++)
		{
			mpara[c][p] = base_para[p]; //VL this copies name, lb, ub, value, fixed, type, delay, display
			mpara[c][p].ival = 0.0;
			mpara[c][p].bayes_mean = 0.0;
			mpara[c][p].bcred_l68 = 0.0;
			mpara[c][p].bcred_l95 = 0.0;
			mpara[c][p].bcred_u68 = 0.0;
			mpara[c][p].bcred_u95 = 0.0;
			mpara[c][p].bcred_median = 0.0;
			mpara[c][p].MLE = 0.0;
			mpara[c][p].MLEl68 = 0.0;
			mpara[c][p].MLEl95 = 0.0;
			mpara[c][p].MLEu68 = 0.0;
			mpara[c][p].MLEu95 = 0.0;
			mpara[c][p].alt = 0;
			mpara[c][p].decis = 0;
			mpara[c][p].delta = 0.0;
			mpara[c][p].pold = 0.0;
			mpara[c][p].padd = 0.0;
			mpara[c][p].altt = 0;
			mpara[c][p].runalt = 0;
			mpara[c][p].runacc = 0;
		}
	}	

	
}
