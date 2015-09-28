//
// implementation of "Fast and Accurate Computation of Binomial Probabilities" by Catherine Loader (July 9, 2000)
//


#include <vector>


namespace filzbach
{
	class sfe : public std::vector<double>
	{
	public:
		sfe();
	};

	double dpois(int x, double lb);
	double dbinom(int x, int n, double p);
}