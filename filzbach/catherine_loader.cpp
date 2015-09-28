#include "stdafx.h"
#include "catherine_loader.h"
#include <assert.h>

namespace filzbach {

	const double PI2=6.283185307179586476925286;

	sfe::sfe()
		: std::vector<double>(16)
	{
		const double log2pi=0.5*log(PI2);
		at(0) = 0.0;
		double lognf = 0.0;
		for (int n=1; n<16; n++) {
			double logn = log((double)n);
			lognf += logn;
			at(n) = lognf+n-log2pi-(0.5+n)*logn;
		}
	}

	double stirlerr(int n)
	{
		static sfe vec_sfe = sfe();
		const double S0 = 1.0/12.0;
		const double S1 = 1.0/360.0;
		const double S2 = 1.0/1260.0;
		const double S3 = 1.0/1680.0;
		const double S4 = 1.0/1188.0;
		if (n<16) return vec_sfe[n];
		double n1 = 1.0/n;
		double n2 = n1*n1;
		if (n>500) return ((S0-S1*n2)*n1);
		if (n>80) return ((S0-(S1-S2*n2)*n2)*n1);
		if (n>35) return ((S0-(S1-(S2-S3*n2)*n2)*n2)*n1);
		return ((S0-(S1-(S2-(S3-S4*n2)*n2)*n2)*n2)*n1);
	}

	double bd0(int x, double np)
	{
		if (fabs(x-np) < 0.1*(x+np)) {
			double v = (x-np)/(x+np);
			double s = (x-np)*v;
			double ej = 2*x*v;
			for (int j=1; ;j++) {
				ej *= v*v;
				double s1 = s + ej/(2*j+1);
				if (s1==s) return s1;
				s = s1;
			}
		}
		return (x*log(x/np)+np-x);
	}

	double dbinom(int x, int n, double p)
	{
		assert((p>=0) && (p<=1));
		assert(n>=0);
		assert((x>=0) && (x<=n));
		if (p==0.0) return x==0 ? 1.0 : 0.0;
		if (p==1.0) return x==n ? 1.0 : 0.0;
		if (x==0) return exp(n*log(1-p));
		if (x==n) return exp(n*log(p));
		double lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x, n*p) - bd0(n-x, n*(1-p));
		return exp(lc)*sqrt(n/(PI2*x*(n-x)));
	}

	double dpois(int x, double lb)
	{
		if (lb==0.0) return x==0 ? 1.0 : 0.0;
		if (x==0) return exp(-lb);
		return exp(-stirlerr(x)-bd0(x,lb))/sqrt(PI2*x);
	}
}