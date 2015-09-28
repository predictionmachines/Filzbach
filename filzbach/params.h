
/* Maximum length of a paramter name */
#define PARAM_NAME_LENGTH 256


// paramter strcut
typedef struct
{  	
	double value;
	double ival;
	double lb;
	double ub;
	char name[PARAM_NAME_LENGTH]; //VL better to use dynamic name allocation: consider (many params)*(many chains)*PARAM_NAME_LENGTH
	int type;			// 0 = additive space, 1 = multiplicative space
	double bayes_mean;
	double bcred_l68;
	double bcred_l95;
	double bcred_u68;
	double bcred_u95;
	double bcred_median;
	double MLE;
	double MLEl68;
	double MLEl95;
	double MLEu68;
	double MLEu95;
	int fixed;
	int alt;
	int decis;
	double delta;
	double pold;
	double padd;
	int altt;			// number of times improved
	int delay;
	int display;
	int runalt;
	int runacc;
	int prioryn;
	double priormean;
	double priorsdev;	
	double rootRhat; // convergence stat for this param (see Gilks et al., 'Markov Chain Monte Carlo in Practice', 1996, page 137)
} param;

extern int paramcount;
extern param* current_para(int currentchain);
extern param* current_para();
extern char* pprintname(int index);
extern void init_chains(int chaincount);
