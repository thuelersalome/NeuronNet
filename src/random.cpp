#include "random.h"

RandomNumbers::RandomNumbers(unsigned long int s)
:seed(s){
     if (seed == 0) {
		 std::random_device rd;
		 seed = rd();
		 }
		 rng = std::mt19937(seed);
	 }
	 
double RandomNumbers::uniform_double(double lower, double upper)
{
	std::uniform_real_distribution<> uniform (lower, upper);
	return uniform(rng);
}
void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper)
{
	for (auto I = vec.begin(); I != vec.end(); I++) {
		*I = uniform_double(lower,upper);
	}
	return;
}
 
double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<> norm (mean, sd);
	return norm(rng);
}
void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd)
{
	for (auto I = vec.begin(); I != vec.end(); I++) {
		*I = normal(mean,sd);
	}
	return;
}

int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<> poi (mean);
	return poi(rng);
}  

void RandomNumbers::poisson(std::vector<int>& vec, double mean)
{ 
   for (auto I = vec.begin(); I != vec.end(); I++) {
		*I = poisson(mean);
	}
	return;
}
