#include "stdafx.h"
#define M_PI       3.14159265358979323846 // On some systems it is also in <math.h>
#include "dwp_stats.h"
#include <float.h>

#include "catherine_loader.h"

/* global variables */

double ncon;
double narray[10003];

///*
// * Lognormal distribution.
// * Based on the formular f(x) = exp(-(ln(x)-mu)^2 / (2*sigma^2)) / (x*(2pi)^0.5) for x > 0
// *
// * x - with x > 0
// * mean - mean value
// * stdev - standard deviation
// */
//double lognormal_density(double x, double mean, double stdev)
//{
//	if (x <= 0)
//	{
//		printf("Error: Paramter x for lognormal_distribution(double x, double mean, double stdev) must be greater than 0.");
//		return nan("");
//	}
//	
//	double lnxm = log(x) - mean;
//	double eln = exp(-(lnxm * lnxm) / (2*pow(stdev, 2)));
//	double xpi = x * pow(2*M_PI, 2) * stdev;
//	double retval = eln / xpi;
//	return retval;
//}

/*
 * Lognormal distribution in logspace.
 */
double lognormal_density(double x, double mode, double stdev)
{
	return normal_density(log(x), log(mode), stdev);
}

///*
// * Returns a lognormal distributed random variable.
// *
// * mean - the mean for the distribution to draw from
// * stdev - the standard deviation of the distribution to draw from in normal space
// */
//double lognormal_draw(double mean, double stdev)
//{
//	return exp(normal_draw(mean, stdev) * stdev + mean);
//}

/*
 * Return a lognormal distributed random variable.
 *
 * mode - the maximum likelihood for the distribution
 * stdev - the standard deviation of the distribution to draw from in log space
 */ 
double lognormal_draw(double mode, double stdev)
{
	return exp(normal_draw(log(mode), stdev));
}

/* 
 * Eponential distribution, describing events occur randonmly over time. 
 *
 * a(lpha) - a > 0, rate paramter
 * x - within 0 <= x < oo
 */
 double exponential_rate(double x, double a)
 {
	if (a <= 0)
	{
		printf("Error: Paramter a for exponential_density(double x, double a) must be greater than 0.");
		return nan("");
	}

	if (x < 0)
	{
		printf("Error: Paramter x for exponential_density(double x, double a) must be greater than or equal 0.");
		return nan("");
	}

	return a * exp(-a * x);
 }

/*
 * Exponential density, using mean instead of the distribution paramter alpha with
 * mean = 1 / alpha.
 */
 double exponential_density(double x, double mean)
 {
	 double a = 1 / mean;
	 return exponential_rate(x, a);
 }

/*
 * Draws a andom number based on the exponential distribution.
 * 
 * a - the rate paramter of the exponential distribution, a > 0
 *
 * TODO: The algorithm is based on the C# build-in uniformyl disributed random number generator.
 *       The code is valid if exp(a) ~ - 1/a + ln()u(0,1), with the unit recangual variante u(0,1);
 *       Need to test if this is valid for rand48 as well, otherwise we have to swith to a different 
 */
 double exponential_rate_draw(double a)
 {
	if (a <= 0)
	{
		printf("Error: Paramter a for exponential_draw(double a) must be greater than 0.");
		return nan("");
	}

	return -log(drand()) / a;
 }
 
/*
 * Draws a andom number based on the exponential distribution using the mean 
 * for the exponential distribution to draw from.
 */
 double exponential_draw(double mean)
 {
	 double a = 1 / mean;
	 return exponential_rate_draw(a);
 }

/*
 * Beta distribution is a continuous probability function distribution 
 * being nonzero only on the interval [0, 1].
 *
 * x - support for x in (0, 1)
 * a(lpha) - shape paramater
 * b(eta) - shape paramter
 */
double beta_density_ab(double x, double a, double b)
{
	// DEAD FUNCTION DO NOT USE
	double bval = beta(a, b);
	return pow(x, a - 1.0) * pow(1.0 -x, b - 1.0) / bval;
}

double beta_density_alpha_beta(double x, double alpha, double beta)
{
	double value;

	// somethings to help, but need more work here
	if(alpha<1.0 && x<0.0010)
		x = 0.0010;

	if(alpha<1.0 && x>0.9990)
		x = 0.9990;

	// follows wikipedia page n pdf of beta distribution

	value = gamma_function(alpha+beta)/(gamma_function(alpha)*gamma_function(beta));

	value *= pow(x,alpha-1.0);

	value *= pow(1.0-x,beta-1.0);

	return value;
}

double beta_density_mean_fvar(double x, double mean, double fvar)
{
	double vcrit = mean*(1.0-mean); // this is the maximum possible variance
	double variance = vcrit*fvar; // actual variance is fvar, times the maximum

	double alpha = mean * (((mean*(1.0-mean))/variance)-1.0);
	double beta = (1.0-mean) * (((mean*(1.0-mean))/variance)-1.0);

	if(alpha<0.0)
	{
		printf("\n error: beta_density_mean_fvar, implied alpha is less than zero \n");
		system("pause");
	}
	if(beta<0.0)
	{
		printf("\n error: beta_density_mean_fvar, implied beta is less than zero \n");
		system("pause");
	}

	return beta_density_alpha_beta(x, alpha, beta);
}

double beta_draw_mean_fvar(double mean, double fvar)
{
	double vcrit = mean*(1.0-mean); // this is the maximum possible variance
	double variance = vcrit*fvar; // actual variance is fvar, times the maximum

	double alpha = mean * (((mean*(1.0-mean))/variance)-1.0);
	double beta = (1.0-mean) * (((mean*(1.0-mean))/variance)-1.0);

	if(alpha<0.0)
	{
		printf("\n error: beta_draw_mean_fvar, implied alpha is less than zero \n");
		system("pause");
	}
	if(beta<0.0)
	{
		printf("\n error: beta_draw_mean_fvar, implied beta is less than zero \n");
		system("pause");
	}

	return beta_draw_alpha_beta(alpha, beta);
}



double beta_draw_alpha_beta(double alpha, double beta)
{
	double maxval=0.0,w,tval;

	for(int ii=0 ; ii<=100 ; ii++)
	{
		double x = 0.0+(double)ii*0.010;
		double val = beta_density_alpha_beta(x, alpha, beta);
		if(val>maxval)
			maxval = val;
	}

	w = 1.0/(10.0*maxval);

	int stop=0;
	while(stop<1)
	{
		tval = random(0.0,1.0);
		double vv = beta_density_alpha_beta(tval, alpha, beta);
		if(random(0.0,1.0)<w*vv)
		{
			stop=1;
		}
	}

	return tval;
}

double beta_draw_mean_alpha(double mean, double alpha)
{
	double beta = alpha*((1.0/mean)-1.0);

	return beta_draw_alpha_beta(alpha,beta);
}



double beta_density_mean_alpha(double x, double mean, double alpha)
{
	double beta = alpha*((1.0/mean)-1.0);

	// follows wikipedia page n pdf of beta distribution

	double value = beta_density_alpha_beta(x, alpha, beta);

	return value;
}




/*
 * Beta distribution using mean and standard deviation. 
 * For calculation of the paramters a and b, mean = a / (a + b) and 
 * var = a*b / ((a+b)(a+b)*(a+b+1)) are used.
 */
double beta_density(double x, double mean, double stdev)
{
	double var = stdev * stdev;
	double msqr = mean * mean;

	double a = (msqr / var) - (1 / var) - (1 / mean);
	double b = (a / mean) - a;

	return beta_density_alpha_beta(x, a, b);
}

/*
 * Draws a random number using the beta distribution.
 *
 * a - used to draw the first gamma random number
 * b - sued to draw the second gamma random number
 */
double beta_draw_ab(int a, int b)
{
	double gamma1 = gamma_draw(a, 1.0);
	double gamma2 = gamma_draw(b, 1.0);
	return gamma1 / (gamma1 + gamma2);
}

/*
 * Probability for number of successes in a sequence
 * of n independed yes or no experiment, each of which yields with 
 * probability p.
 *
 * x - the number of draws to get the probability for.
 * n - number of trials, n >= 0
 * p - success probability in each trial, p in [0,1]
 */
double binomial_density(int x, int n, double p)
{
	CHECK (n < 0, "Error: Paramter n for binamial_density(int x, int n, double p) must be greater or equal than 0.");
	CHECK (x < 0 || x > n, "Error: Paramter x for binamial_density(int x, int n, double p) must be within [0, n].");
	CHECK (p < 0.0 || p > 1.0, "Error: Paramter p for binamial_density(int x, int n, double p) must be within [0, 1].");
	return filzbach::dbinom(x,n,p);
}

/*
 * Draws a random number out of the binamial distribution.
 *
 * n - number of trials, n >= 0
 * p - success probability in each trial, p in [0,1]
 */
double binomial_draw(int n, double p)
{
	if (n < 0)
	{
		printf("Error: Paramter n for binamial_density(int x, int n, double p) must be greater or equal than 0.");
		return nan("");
	}

	if (p < 0.0 || p > 1.0)
	{
		printf("Error: Paramter p for binamial_density(int x, int n, double p) must be within [0, 1].");
		return nan("");
	}	

	double res = 0.0;

	for (int i = 0; i < n; i++)
	{
		if (drand() < p)
			res++;
	}

	return res;
}

/*
 * Negative binomial density.
 *
 * k - 
 * r - r > 0 -  number of failures until the experiment is stopped 
 * p - in (0,1) - success probability in each experiment (real)
 */
double negative_binomial_density_p(int k, int r, double p)
{
	if (r <= 0.0)
	{	
		printf("Error: Paramter n for negative_binomial_draw(double r, double p) must be greater than 0.");
		return nan("");
	}

	if (p < 0.0 || p > 1.0)
	{
		printf("Error: Paramter p for negative_binomial_draw(int r, double p) must be within [0, 1].");
		return nan("");
	}

	int a = k + r - 1;
	int b = r - 1;

	return factorial(a) / (factorial(a-b) * factorial(b));
}

/*
 * Negative binomial density using mean and standard deviation.
 *
 * Calling the common negative binomial distribution with paramters k, r, and p
 * calculating the required paramters based on the forumlars r * ( p / (1-p) ) for mean
 * and r * (p / (1 - p)(1 - p)) for the variance.
 * 
 */
double negative_binomial_density(int k, double mean, double stdev)
{
	double p = 1 - (mean / (stdev*stdev));
	double r = ((mean / p) * (1 - p));
	return negative_binomial_density_shape(k, r, p);
}

/*
 * Negative binomial density using mean and standard deviation.
 *
 * Calling the common negative binomial distribution with paramters k, r, and p
 * calculating the required paramters based on the forumlars r * ( p / (1-p) ) for mean
 * and r * (p / (1 - p)(1 - p)) for the variance.
 * 
 */
double negative_binomial_density_stdev(int k, double mean, double stdev)
{
	double p = 1 - (mean / (stdev*stdev));
	double r = ((mean / p) * (1 - p));
	return negative_binomial_density_shape(k, r, p);
}


/*
 * Negative binomial density.
 *
 * ss - number of obersvation
 * k shape 
 *
 * TODO: Paramter'code' has been temporarily disabled
 */
double negative_binomial_density_shape(int ss, double r, double shape) //, int code( 
{
	int code = 0; 
	static int nbsetup=0;
	static double (*negbinom_pzero)[10002];

	/* iterative procedure lifted from ecological detective */
	int xx;
	double prob, prob_prev;
	double s, tk, tmean;
	int kk, mm;

	if (nbsetup<1 && code>0)
	{
		printf("\n setting up storage array for negative binomial \n");
		negbinom_pzero = new double[10002][10002];
		/* set up array for power part */
		for(kk=1 ; kk<=10000 ; kk++){
			for(mm=1 ; mm<=10000 ; mm++)
			{
				tk = (3.0 / 10000.0) * (double)(kk+1);
				tmean = (5.0 / 10000.0) * (double)(mm+1);
				/* p zero  */
				negbinom_pzero[kk][mm] = pow(tk/(tk+tmean),tk);
			}
		}
		nbsetup = 1;
	}

	if (code>0 && r <= 4.95 && r>=0.00060 && shape<=2.95 && shape>=0.00040){
		kk = (int)((shape/3.0)*10000.0);
		mm = (int)((r/5.0)*10000.0);
		// printf("\n negbinom array r %f mm %d k %f kk %d",r,mm,k,kk);
		prob = negbinom_pzero[kk][mm];
	}
	else
	{
		prob = pow(shape/(shape+r),shape);
	}

	if(ss>0)
	{
		for(xx=1 ; xx<=ss ; xx++){
			s=(double)xx;
			prob_prev = prob;
			prob = ((s+shape-1.0)/s)*(r/(shape+r))*prob_prev;
		}
	}

	return prob;
}

/*
 * Draws a random number out of the negative binamial distribution.
 * The form of the distribution is not the number of trials to the n-th success. 
 *
 * n - number of trials, n >= 0
 * p - success probability, p in [0,1]
 *
 * TODO: the description of calculation y in "y = y * p / (1.0 - p)" was rarther vague, this reuqires some verification for correctness.
 */
double negative_binomial_draw_p(double n, double p)
{
	if (n < 0.0)
	{	
		printf("Error: Paramter n for negative_binomial_draw(int n, double p) must be greater than 0.");
		return nan("");
	}

	if (p < 0.0 || p > 1.0)
	{
		printf("Error: Paramter p for negative_binomial_draw(int n, double p) must be within [0, 1].");
		return nan("");
	}
		

	// draw Gamma random variable to determine Y.
	double y = gamma_draw(n, 1.0);
	y = y * p / (1.0 - p); 

	// draw Poisson random variable with mean of Y
	return poisson_draw(y);
}

/*
 * Draws a random number out of the negative binamial distribution.
 * Uses a given mean and standard deviation.
 */
double negative_binomial_draw(double mean, double stdev)
{
	double p = 1 - (mean / (stdev*stdev));
	double r = ((mean / p) * (1 - p));
	return negative_binomial_draw_p(r, p);
}

/*
 *
 */
void testing_neg_binom()
{
	int kk,mm,ss;
	double mean,k;

	for (mm=1 ; mm<=10 ; mm+=2)
	{
		for (kk=1 ; kk<=10 ; kk+=2)
		{
			for(ss=0 ; ss<=6 ; ss+=2)
			{
				mean = 0.10 * (double)mm;
				k = 0.50 * (double)kk;
				printf("\n mean %f k %f ss %d",mean,k,ss);
				printf("\n act %f from_array %f",negative_binomial_density_shape(ss,mean,k/*,0*/),negative_binomial_density_shape(ss,mean,k/*,1*/));
			}
		}
	}
}

/*
 *
 */
double normal_density(double n, double mean, double stdev)
{
	static int nsetup=0;

	double temp,ratio,tval,tmean,tstdev,tvar;
	int rr;
	double mult;

	if(nsetup<1)
	{
		ncon=sqrt(2.0*3.1415927);

		/* set up array that stores density vs ratio of deviation to stdev */
		for(rr=0 ; rr<=10000 ; rr++)
		{
			ratio=(((double)rr)/10000.0)*12.0; /* only go up to 60 stdevs -- that's a lot!! */
			tmean=100.0;
			tstdev=100.0;
			tval=tmean+(ratio*tstdev);
			tvar=tstdev*tstdev;
			// temp=1.0/sqrt(2.0*3.1415927*tvar); /* constant not independent of mean and val */
			temp=exp(((-1.0)*pow(tval-tmean,2.0))/(2.0*tvar)); /* ...but this part is (i think!) */
			narray[rr]=temp;
			// if(rr%10==0) printf("\n %d --> %f ",rr,temp);
		}
	}

	nsetup=1;

	/* calc ratio of (val-mean) to stdev */
	ratio=(fabs(n-mean))/stdev;
	ratio=(ratio/12.0)*10000.0;
	temp = 0;
	if ((ratio<=9999) && (ratio>=0))
	{
		rr=(int)(ratio);

		mult=narray[rr];
		temp=1.0/(ncon*stdev);
		// temp = 1.0/(sqrt(2.0*3.1415927)*stdev); /* constant not independent of mean and val */
		temp*=mult;
	}
	/*
	printf(" | q:%f",temp);
	mult = exp(((-1.0)*pow(val-mean,2.0))/(2.0*stdev*stdev));
	temp = 1.0/(sqrt(2.0*3.1415927)*stdev); 
	temp*=mult;
	printf(",s:%f |",temp);
	*/

	if(temp<0.00000000010)
		temp=0.00000000010;

	return temp;
}

double gaussian(double x, double mean, double sigma)
{
	static int gsetup=0;
	static double height[201];

	double temp;
	double edev1,edev2;
	int i,i1,i2;

	if(gsetup<1)
	{
		// do 100 different values of scaled deviation from mean and store in static array for lookup hereafter
		for(i = 0 ; i<=200 ; i++)
		{
			double scdev = 10.0 * (((double)i)/200.0);
			height[i] = exp((-1.0) * pow(scdev,2.0));
		}
		gsetup=1;
	}

	// calculate scdev for this lookup
	double scdev = (x-mean)/sigma;
	if(scdev<0.0)
		scdev*=-1.0;

	// now lookup if scdev<5.0 or calculate raw otherwise
	if(scdev<10.0)
	{
		i1 = (int)((scdev/10.0)*200.0);
		edev1 = ((double)i1)*(10.0/200.0);
		i2 = i1+1;
		edev2 = ((double)i2)*(10.0/200.0);
		// interpolate
		temp = ((scdev-edev1)/(edev2-edev1))*height[i2] + ((edev2-scdev)/(edev2-edev1))*height[i1];
	}
	else
	{
		temp = exp((-1.0) * pow(scdev,2.0));
	}

	return temp;
}

/*
 * Box-Muller transform
 */
bool normal_draw_phase2 = false;
double normal_draw(double mean, double stdev)
{
	static double V2, fac;

	double S, Z, U1, U2, V1;

	if (normal_draw_phase2)
	{
		Z = V2 * fac;
		normal_draw_phase2 = false;
	}
	else
	{
		do {
			U1 = drand();
			U2 = drand();
			V1 = (2.0*U1)-1.0;
			V2 = (2.0*U2)-1.0;
			S = (V1*V1)+(V2*V2);
		} while (S>=1.0);

		fac = sqrt(-2.0 * (log(S)/S));
		Z = V1 * fac;
		normal_draw_phase2 = true;
	}

	Z*=stdev;
	Z+=mean;

	return Z;
} 

/*
 * Gamma density.
 *
 * n observations
 * mean
 * stdev standard deviation
 */
extern double gamma_density(double n, double mean, double stdev)
{
	double q = mean / stdev;
	double shape = q * q;

	return gamma_density_shape(n, mean, shape);
}

/*
 * Gamma density
 *
 * n - observations
 * shape - shape paramter
 * sclae - scale paramter
 */
double gamma_density_scale(double n, double shape, double scale)
{
	double mean = shape * scale;
	return gamma_density_shape(n, mean, shape);
}

/*
 * Gamma density.
 *
 * n observations
 * mean
 * shape 
 */

// mean -> mean
// z = observation
// n shape paramte
double gamma_density_shape(double n, double mean, double shape)
{
	double cap_rho;
	double ck[20];
	int k;
	double cknk_sum;
	double f_z;
	double a;
	int numit;
	double comp,comp2;
	int ii;

	/* Christian -- this little section looks a bit worrying but I can't remember whether it really matters or not! */

	/* correction because when n gets very small this approx is not good */
	if(shape<0.0020)
		shape=0.0020;
	/* same for z value */
	if(n<0.0020)
		n=0.0020;

	/* this is the trick with the mean */

	a = shape / mean;

	/* write ck values for use a little later*/

	ck[1]=1.0;
	ck[2]=0.5772156649015329;
	ck[3]=-0.6558780715202538;
	ck[4]=-0.0420026350340952;
	ck[5]=0.1665386113822915;
	ck[6]=-0.0421977345555443;
	ck[7]=-0.0096219715278770;
	ck[8]=0.0072189432466630;
	ck[9]=-0.0011651675918591;
	ck[10]=-0.0002152416741149;
	ck[11]=0.0001280502823882;
	ck[12]=-0.0000201348547807;
	ck[13]=-0.0000012504934821;
	ck[14]=0.0000011330272320;
	ck[15]=-0.0000002056338417;
	ck[16]=0.0000000061160950;
	ck[17]=0.0000000050020075;
	ck[18]=-0.0000000011812746;
	ck[19]=0.0000000001043427;

	/* calc constant cap_rho(n) */

	/* recurrence relation thing (see hilborn and mangel) */
	numit = (int)shape;
	comp = (double)numit;
	comp2 = shape-comp;
	//so, comp2 is now in the range 0-1, which when addded to comp, gives our original n (here, n is notated as shape, sorry)
	/* if n is a whole number then potential prob -- need to correct as follows */
	if (fabs(comp-shape)<0.001)
	{
		comp2 += 1.0;
		numit -= 1;
	}


	/* calc fraction part */
	cknk_sum=0.0;
	for (k=1 ; k<=19 ; k++)
	{
		cknk_sum += ck[k] * pow(comp2,(double)k);
	}

	cap_rho = 1.0 / cknk_sum;



	for (ii=1 ; ii<=numit ; ii++)
	{
		cap_rho *= shape-(1.0*((double)ii));
	}

	/* calc f_z */

	// notation not great here. Matching to Hilborn and Mangel (page 77) we have their f_z = our f_z, their cap_rho_n = our cap_rho, their n = our shape, their a = our a, their z = our n.
	// so, above, cap_rho of n means 
	f_z = (pow(a,shape)/(cap_rho)) * exp((-1.0)*a*n) * pow(n,shape-1);

	return f_z;
}

double gamma_draw(double mean, double stdev)
{
	return nan("");
}

/*
 * Draws a single random number out of the gamma distribution with 
 * a given scale and shape.
 */
double gamma_draw_shape(int scale, double shape)
{
	double res = 0.0;

	for (int i = 0; i < scale; i++)	
		res += -log(drand()) / shape;	

	return res;
}

/*
 * Draws a single random number out of the gamma distribution with 
 * a given scale and mean.
 */
double gamma_draw_mean(int scale, double mean)
{
	double shape = mean / scale;
	return gamma_draw(scale, shape);
}

/*
 * Draws a single random number out of the gamma distribution with 
 * a given scale and standard deviation.
 */
extern double gamma_draw_stdev(int scale, double stdev)
{
	double shape = (stdev * stdev)  / (scale * scale);
	return gamma_draw(scale, shape);
}

double gamma_function(double n)
{
	
	/* write ck values for use a little later*/
	double ck[21];

	ck[1]=1.0;
	ck[2]=0.5772156649015329;
	ck[3]=-0.6558780715202538;
	ck[4]=-0.0420026350340952;
	ck[5]=0.1665386113822915;
	ck[6]=-0.0421977345555443;
	ck[7]=-0.0096219715278770;
	ck[8]=0.0072189432466630;
	ck[9]=-0.0011651675918591;
	ck[10]=-0.0002152416741149;
	ck[11]=0.0001280502823882;
	ck[12]=-0.0000201348547807;
	ck[13]=-0.0000012504934821;
	ck[14]=0.0000011330272320;
	ck[15]=-0.0000002056338417;
	ck[16]=0.0000000061160950;
	ck[17]=0.0000000050020075;
	ck[18]=-0.0000000011812746;
	ck[19]=0.0000000001043427;

	/* recurrence relation thing (see hilborn and mangel) */
	int numit = (int)n;
	double n_I = (double)numit;
	double n_f = n-n_I;
	//so, n_f is now in the range 0-1, which when addded to n_I, gives our original n 
	/* if n is already whole number then potential prob -- need to correct as follows */
	if (fabs(n_I-n)<0.001)
	{
		n_f += 1.0;
		numit -= 1;
	}

	/* calc rho(n_f) */
	double cknk_sum=0.0;
	for (int k=1 ; k<=19 ; k++)
	{
		cknk_sum += ck[k] * pow(n_f,(double)k);
	}
	double cap_rho = 1.0 / cknk_sum;

	// then do the easy multiplication bit
	for (int ii=1 ; ii<=numit ; ii++)
	{
		cap_rho *= n-(1.0*((double)ii));
	}

	return cap_rho;
	
}




/*
 *
 */
double factorial(int k)
{
	static int ldone=0;
	static double faclist[201];

	double val=1;
	int ii,vv;
	
	if (ldone<1)
	{
		ldone=1;
		faclist[0]=1.0;
		for (vv=1 ; vv<=200 ; vv++)
		{
			val=1;
			for( ii=2 ; ii<=vv ; ii++)
			{
				val*=(double)ii;
			}
			faclist[vv]=val;
		}
	}

	if(k<=200)
	{
		return faclist[k];
	}
	else
	{
		val=1;
		for(ii=2 ; ii<=k ; ii++)
		{
			val*=(double)ii;
		}
	}

	return val;
}

/*
 * Poisson density.
 *
 * n - number of observations.
 * lambda - excptected number of occurences.
 */
double poisson_density(int n, double lambda)
{
	return filzbach::dpois(n,lambda);
}


/*
 * lambda - excptected number of occurences.
 */
int poisson_draw(double lambda)
{
	double ell = exp((-1.0)*lambda);
	int k=0;
	double p=1.0;

	if (lambda<30.0)
	{
		while(p >= ell)
		{
			k++;
			p *= drand();
		}
		return k-1;
	}
	else
	{
		return (int)normal_draw(lambda,sqrt(lambda));
	}
}


/*
 * 
 */
double spearman_rank(double xx[], double yy[], int n)
{
	double *xrank = new double[n]; 
	double *yrank = new double[n];
	short int *used_x = new short int[n];
	short int *used_y = new short int[n];
	
	for (int p = 0; p < n; p++) 
	{
		used_x[p] = 0;
		used_y[p] = 0;
	}
	
	/* ranks for xvals */
	for(int x1 = 0; x1 < n; x1++)
	{
		int pos = -1;
		double min = 999999.0;

		for(int x2 = 0; x2 < n; x2++)
		{			
			if(used_x[x2] < 1 && xx[x2] < min)
			{
				min = xx[x2];
				pos = x2;
			}
		}

		if (pos > -1)
		{
			xrank[pos] = x1 + 1;
			used_x[pos] = 1;
		}
	}		

	/* ranks for yvals */
	for(int y1 = 0; y1 < n; y1++)
	{
		int pos = -1;
		double min = 999999.0;

		for(int y2 = 0; y2<n; y2++)
		{
			if(used_y[y2] < 1 && yy[y2] < min)
			{
				min = yy[y2];
				pos = y2;  
			}
		}

		if (pos > -1)
		{
			yrank[pos] = y1 + 1;
			used_y[pos] = 1;
		}
	}

	double sumdifs=0.0;
	for (int i=0 ; i<n ; i++)
	{
		double dif = (double)xrank[i] - (double)yrank[i];
		sumdifs += dif*dif;
	}

	return 1.0 - ((6.0 * sumdifs) / (pow((double)n, 3.0) - (double)n));
}

/*
 *
 */
void test_spearman()
{
	double xvals[6]={0.50,1.0,1.50,2.0,2.50,3.0};
	double yvals[6]={3.0,2.50,2.0,1.50,1.0,0.50};
	double rs;

	printf("\n testing SPEARMAN \n");
	rs = spearman_rank(xvals,yvals,6);
	printf(" rs = %f \n",rs);
} 

/*
 *
 */
double logistic(double emm)
{
	static int logit_setup=0;
	static double *logit_array;
	int bin;
	double etm,rval;

	if(logit_setup<1)
	{
		/* set up array */
		logit_array = new double[20002];

		for(bin = 0 ; bin<=20000 ; bin++)
		{
			etm = (double)bin - 10000.0;
			etm /= 1000.0;
			logit_array[bin]=1.0/(1.0+exp((-1.0)*etm));
		}
		logit_setup=1;
	}

	bin = (int)((emm*1000.0)+10000.0);

	if (bin>0 && bin<20000)
	{
		rval = logit_array[bin];
	}
	else
	{
		rval = 1.0 / (1.0+exp((-1.0)*emm));
		// printf("\t logit outside bin range ");
	}

	return rval;
}

/*
 *
 */
void testing_logit()
{
	double emm;
	int kk;
	double tval;

	for (kk=1 ; kk<=100 ; kk++)
	{
		emm = -5.0+(0.10*(double)kk);
		tval = 1.0/(1.0+exp((-1.0)*emm));
		printf("\n true logit %f",tval);
		printf("\t logit returned %f",logistic(emm));
	}
}



/*
 * Random floating point number, uniform distribution from [0, 1)
 */
double drand()
{
	double r;
	#pragma omp critical
	{
		r = ((double)rand()) / (1.0 + (double)RAND_MAX);
	}
	return r;
} 

/*
 *	Gamma method based on Lanczos approximation.
 */
double gamma(double x)
{
	// DEAD FUNCTION DO NOT USE
	const int G = 7;

	double coef[9];
	coef[0] = 0.99999999999980993;
	coef[1] = 676.5203681218851;
	coef[2] = -1259.1392167224028;
	coef[3] = 771.32342877765313;
	coef[4] = -176.61502916214059;
	coef[5] = 12.507343278686905;
	coef[6] = -0.13857109526572012;
	coef[7] = 9.9843695780195716e-6;
	coef[8] = 1.5056327351493116e-7;

	if ( x < 0.5)
		return M_PI / (sin(M_PI * x) * gamma(1.0 - x));

	x -= 1.0;

	double y = coef[0];
	for (int i = 1; i < G + 2; i++)
		y += coef[i] / (x + 1.0 * i);

	double z = x + (G + 0.5);

	return sqrt(2 * M_PI) * pow(z, x + 0.5) * exp(-z) * y;
}

/*
 * Beta function, Euler integral of the first kind.
 */
double beta(double x, double y)
{
	// DEAD FUNCTION DO NOT USE
	return gamma(x) * gamma(y) / gamma(x + y);
}

/***********************************************/
double random(double lb, double ub)
{
	return lb + drand()*(ub-lb);
}

/***********************************************/
int random_integer(int lb, int ub)
{
	return lb + (int)( ((double)(ub-lb+1))*drand() );
}

/**********************************************
 * Find 'num'-th largest value in the section
 * value[lid]..value[uid]
*/

double large(double value[], int lid, int uid, int num)
{
	int nn,ii,used[100001],maxid=0,hits;
	double maxval = log(0.0); // -Infinity

	if(uid>100000)
	{
		printf("\n error in 'large' -- can only deal with lists with indices up to 100,000 -- sorry \n");
		system("pause");
		return -1;
	}
	else
	{
		for(ii=lid ; ii<=uid ; ii++)
			used[ii]=0;

		for(nn=1 ; nn<=num ; nn++)
		{
			hits=0;
			for(ii=lid ; ii<=uid ; ii++)
			{
				if(used[ii]<1)
				{
					hits++;
					if(value[ii]>maxval || hits==1)
					{
						maxval = value[ii];
						maxid = ii;
					}
				}
			} // ii loop
			used[maxid]=1;
		} // nn loop
	}

	return value[maxid];

}

/***********************************************/
double large(float value[], int lid, int uid, int num)
{
	int nn,ii,used[100001],maxid=0,hits;
	double maxval = log(0.0); // -Infinity

	if(uid>100000)
	{
		printf("\n error in 'large' -- can only deal with lists with indices up to 100,000 -- sorry \n");
		system("pause");
		return -1;
	}
	else
	{
		for(ii=lid ; ii<=uid ; ii++)
			used[ii]=0;

		for(nn=1 ; nn<=num ; nn++)
		{
			hits=0;
			for(ii=lid ; ii<=uid ; ii++)
			{
				if(used[ii]<1)
				{
					hits++;
					if(value[ii]>maxval || hits==1)
					{
						maxval = value[ii];
						maxid = ii;
					}
				}
			} // ii loop
			used[maxid]=1;
		} // nn loop
	}

	return value[maxid];

}

/***********************************************/
double large(int value[], int lid, int uid, int num)
{
	int nn,ii,used[100001],maxid=0,hits;
	int maxval = value[0]-1;
	for (int i=1; i<lid; i++) if (maxval>=value[i]) maxval=value[i]-1;

	if(uid>100000)
	{
		printf("\n error in 'large' -- can only deal with lists with indices up to 100,000 -- sorry \n");
		system("pause");
		return -1;
	}
	else
	{
		for(ii=lid ; ii<=uid ; ii++)
			used[ii]=0;

		for(nn=1 ; nn<=num ; nn++)
		{
			hits=0;
			for(ii=lid ; ii<=uid ; ii++)
			{
				if(used[ii]<1)
				{
					hits++;
					if(value[ii]>maxval || hits==1)
					{
						maxval = value[ii];
						maxid = ii;
					}
				}
			} // ii loop
			used[maxid]=1;
		} // nn loop
	}

	return value[maxid];

}

void initialize_stat()
{
	normal_draw_phase2=false;
}
