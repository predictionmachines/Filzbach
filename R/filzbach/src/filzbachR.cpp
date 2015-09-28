// Interop.cpp : Defines the exported functions for the DLL application.
//

#include "filzbach.h"
#include "R.h"
#include "Rinternals.h"
//#include "filzbachR.h"

#include <functional>


void R_likelihood(
	SEXP call_likelihood, // a LANGSXP pairlist prepared to call R log-likelihood function
	long sample_size, // a number to pass to set_metr_number_ok()
	SEXP param_defs, // a list describing filzbach parameters
	SEXP env) // an R environment to call R log-likelihood in
{
	SEXP names = getAttrib(param_defs, R_NamesSymbol);
	int c_params = length(names);
	if (!isString(names) || c_params<=0) error("filzbachRunMCMC: (*) no param names");
	SEXP t = CDR(call_likelihood);
	for (int i=0; i<c_params; i++)
	{
		const char* name = CHAR(STRING_ELT(names,i));
		SEXP def = VECTOR_ELT(param_defs,i);
		if (length(def)==6)
			REAL(CAR(t))[0] = cv(name);
		else for (int j=0; j<REAL(def)[6]; j++)
			REAL(CAR(t))[j] = cv(name, j);
		t = CDR(t);
	}
	SEXP result = eval(call_likelihood, env);
	if (!isReal(result)) error("filzbachRunMCMC: (*) 'likelihood' result must be a number");
	set_metr_ltotnew(REAL(result)[0]);
	set_metr_number_ok(sample_size);
}

std::function<void()> lambda;
void relay() {lambda();}

// param_defs is a list of parameters

extern "C" SEXP filzbachRunMCMC(SEXP args) // (burnin, eststeps, likelihood, sample_size, param_defs, thinning, environment) -> (matrix of the resulting bayes chain)
{
	SEXP ans, dims, dimnames;
	if (!isList(args)) error("filzbachRunMCMC: argument not a List. You sould call this function using .External");
	if (!isReal(CADR(args)) && !isInteger(CADR(args))) error("filzbachRunMCMC: 'burnin' must be a numeric value");
	if (!isNumeric(CADDR(args))&& !isInteger(CADDR(args))) error("filzbachRunMCMC: 'eststeps' must be a numeric value");
	SEXP t = CDR(CDDR(args));
	SEXP likelihood = CAR(t);
	if (!isFunction(likelihood)) error("filzbachRunMCMC: 'likelihood' must be a function");
	if (!isReal(CADR(t)) && !isInteger(CADR(t))) error("filzbachRunMCMC: 'sample_size' must be a number");
	long sample_size = isReal(CADR(t)) ? (long)REAL(CADR(t))[0] : (long)INTEGER(CADR(t))[0];
	SEXP param_defs = CADDR(t);
	if (!isNewList(param_defs)) error("filzbachRunMCMC: 'param_defs' must be a list");
	SEXP param_names = getAttrib(param_defs, R_NamesSymbol);
	int c_params = length(param_names);
	if (!isString(param_names) || c_params<=0) error("filzbachRunMCMC: (*) no param names");
	t = CDDR(t);
	SEXP thinning = CADR(t);
	if (!isNumeric(thinning)) error("filzbachRunMCMC: 'thinning' must be a numeric value");
	SEXP env = CADDR(t);
	if (!isEnvironment(env)) error("filzbachRunMCMC: 'environment' must be an environment");

	initialize_filzbach();
	for (int i=0; i<c_params; i++)
	{
		SEXP def = VECTOR_ELT(param_defs,i);
		const char* name = CHAR(STRING_ELT(param_names,i));
		if (!isReal(def) || length(def)<6 || length(def)>7) error("filzbachRunMCMC: invalid definition for parameter %d (%s)", i, name);
		double* v = REAL(def);
		if (length(def)==6)
			parameter_create(name,v[0],v[1],v[2],(int)v[3],(int)v[4],(int)v[5]);
		else
			parameter_create_vector(name,v[0],v[1],v[2],(int)v[3],(int)v[4],(int)v[5],(int)v[6]);
	}
	// prepare a LANGSXP pairlist to call R log-likelihood function
	SEXP call_likelihood;
    PROTECT(call_likelihood = allocList(c_params+1));
    SET_TYPEOF(call_likelihood, LANGSXP);
    SETCAR(call_likelihood, likelihood);
	t = CDR(call_likelihood);
	for (int i=0; i<c_params; i++)
	{
		SEXP def = VECTOR_ELT(param_defs,i);
		if (length(def)==6)
			SETCAR(t, allocVector(REALSXP,1));
		else
			SETCAR(t, allocVector(REALSXP,(int)REAL(def)[6]));
		t = CDR(t);
	}
	lambda = [call_likelihood,sample_size, param_defs, env](){R_likelihood(call_likelihood, sample_size, param_defs, env);};
	pfn_likelihood=relay;
	set_thinning(asInteger(thinning));
	runmcmc(asInteger(CADR(args)),asInteger(CADDR(args)),100,100);

	UNPROTECT(1);
	// compute resulting matrix dimensions
	int param_len=0;
	for (int i=0; i<c_params; i++)
	{
		SEXP def = VECTOR_ELT(param_defs,i);
		param_len += (length(def)==6) ? 1 : REAL(def)[6]; 
	}
	int chain_len=0;
	while (params_from_bayes_list(0,chain_len)==0) chain_len++;
	// allocate matrix
	PROTECT(ans = allocVector(REALSXP, chain_len*param_len));
	// set dims
	PROTECT(dims = allocVector(INTSXP, 2));
	INTEGER(dims)[0] = chain_len;
	INTEGER(dims)[1] = param_len;
	setAttrib(ans, R_DimSymbol, dims);
	// set dimnames
	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, R_NilValue);
	SET_VECTOR_ELT(dimnames, 1, allocVector(STRSXP, param_len));
	SEXP colnames = VECTOR_ELT(dimnames, 1);
	int col = 0;
	for (int i=0; i<c_params; i++)
	{
		SEXP def = VECTOR_ELT(param_defs,i);
		if (length(def)==6)
			SET_STRING_ELT(colnames, col++, STRING_ELT(param_names,i));
		else for (int j=0; j<REAL(def)[6]; j++)
		{
			char buf[128];
			snprintf(buf,sizeof(buf),"%s[%d]",CHAR(STRING_ELT(param_names,i)),j);
			SET_STRING_ELT(colnames, col++, mkChar(buf));
		}
	}
	setAttrib(ans,R_DimNamesSymbol, dimnames);
	// set data
	for (int j=0; j<chain_len; j++)
	{
		params_from_bayes_list(0,j);
		int col=0;
		for (int i=0; i<c_params; i++)
		{
			SEXP def = VECTOR_ELT(param_defs,i);
			const char* name = CHAR(STRING_ELT(param_names,i));
		if (length(def)==6)
			REAL(ans)[(col++)*chain_len+j] = cv(name);
		else for (int k=0; k<REAL(def)[6]; k++)
			REAL(ans)[(col++)*chain_len+j] = cv(name,k);
		}
	}
	UNPROTECT(3);
	return ans;
}

