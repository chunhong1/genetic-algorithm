#include <iostream>
#include <conio.h>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <string>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <windows.h>
#include <cmath>
#include <limits>
using namespace std;

/******************************************************** CODE GUIDELINE ********************************************************
* Parameter Settings

- Benchmark Function
- External Function
* Selection Function
* Crossover Function
* Mutation Function
* Replacement Function

main()
- GA Automation
	-> Benchmark Automation
		-> Experiment Automation
			-> Generate Population
			-> Fitness Evaluation 1
			-> Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
				*> Selection Operation
				*> Crossover Operation
				*> Mutation Operation
				-> Fitness Evaluation 2
				*> Replacement Operation
******************************************************** CODE GUIDELINE ********************************************************/

//------------------------------------------------------------------------------------------------------------------------------
// Parameter Settings
//------------------------------------------------------------------------------------------------------------------------------
#define MINI_PROJECT -2 							// Project -> 1 | Assignment -> 0 | Demo -> -1 | Manual -> -2

//////////////////////////
/* GA_COMBINATION Sheet */
//////////////////////////
/* Selection	- [0] */
// Best 			[0][1] Manual 			[0][1]

// Roulette Wheel	[0][0]
// Tournament		[0][1]
// Linear Ranking	[0][2]
// XXX				[0][3]
//////////////////////////
/* Crossover	- [1] */
// Best 			[1][1] Manual 			[0][1]

// Uniform			[1][0]
// Shuffle			[1][1]
// Arithmetic		[1][2]
// XXX				[1][3]
//////////////////////////
/* Mutation		- [2] */
// Best 			[2][1] Manual 			[0][1]

// Reversing		[2][0]
// Random			[2][1]
// Simple Inversion	[2][2]
// XXX				[2][3]
//////////////////////////
/* Replacement 	- [3] */
// Best 			[3][1] Manual 			[0][1]

// Weak Parent		[3][0]
// Both Parent		[3][1]
// Tournament		[3][2]
// XXX				[3][3]
//////////////////////////

const double dcp = 0.7, dmp = 0.01;				// Crossover Probability, Mutation Probability
const int tournamentSize = 5;					// Tournament Selection Size
int dynamicTournamentSize = tournamentSize;
const float highDiverseThreshold = 0.04;		// threshold for high diversity
const float lowDiverseThreshold = 0.01;			// threshold for low diversity
//------------------------------------------------------------------------------------------------------------------------------
// #if MINI_PROJECT == 1
// 	#define TECHNIQUE 4							// No. Of Techniques
// 	string GA = "0000000000000000";				// GA Binary Combinations
// #elif MINI_PROJECT == 0 || MINI_PROJECT == -1 || MINI_PROJECT == -2
// 	#define TECHNIQUE 2
// 	string GA = "00000000";
// #else
// 	#define EXIT
// 	#define TECHNIQUE 1
// 	string GA = "";
// #endif

#if MINI_PROJECT == 0 || MINI_PROJECT == -1 || MINI_PROJECT == -2
	#define TECHNIQUE 2
	string GA = "00000000";
#else
	#define EXIT
	#define TECHNIQUE 1
	string GA = "";
#endif

#define getrandom(min, max) (static_cast<long long>(rand()) * (max - min + 1) / RAND_MAX) + min
#define gen 2000								// number of iterations (number of generations)
#define pSize 40								// number of chromosomes (population size)
#define dimension 30							// number of bits (dimension size)

double chromosome[pSize][dimension];			// chromosome
double paroff[4][dimension];					// parent and offspring
double fit[pSize];								// fitness value for each chromosome
double r = 0, gcp = 0, gmp = 0;
int crb = 0, mb1 = 0, mb2 = 0;
int rp1 = 0, rp2 = 0;
double mb1v = 0, mb2v = 0;
double fv = 0, sumFit = 0;
double fit1 = 0, fit2 = 0;
double tfit[4];

int lFvIndex = 0;
double lFv = 0;
double lowestGene[dimension];
double lowestGeneFV = 0;

const double pi = 2 * asin(1.0);
int GA_COMBINATION[4][TECHNIQUE];
int BENCHMARK = 1;
string benchmarkFunction = "";
int rangeMin = 0, rangeMax = 0, rangeDiv = 1000;

//------------------------------------------------------------------------------------------------------------------------------
// Benchmark Function
//------------------------------------------------------------------------------------------------------------------------------
// No.1 - Sphere Function +-5.12
double Sphere(double a[])
{
	for (int j = 0; j < dimension; j++)
	{
		fv = pow(a[j], 2);
		sumFit = sumFit + fv;
	}

	return sumFit;
}

// No.2 - Ackley Function +-32.768
// Done By: Yeap Chun Hong 2206352
double Ackley(double a[])
{
	double sum1 = 0, sum2 = 0;

	for (int j = 0; j < dimension; j++)
	{
		sum1 += pow(a[j], 2);
		sum2 += cos(2 * pi * a[j]);
	}

	double term1 = -20 * exp(-0.2 * sqrt(sum1 / dimension));
	double term2 = exp(sum2 / dimension);

	sumFit = term1 - term2 + 20 + exp(1);
	return sumFit;
}

// No.3 - Rastrigin Function +-5.12
// Done By: Yeap Chun Hong 2206352
double Rastrigin(double a[])
{
	for (int j = 0; j < dimension; j++)
	{
		fv = (pow(a[j], 2)) - (10 * cos(2 * pi * a[j]));
		sumFit = sumFit + fv;
	}

	sumFit += 10 * dimension;
	return sumFit;
}

// No.4 - Zakharov Function -5, +10
// Done By: Brandon Ting En Junn 2101751
double Zakharov(double a[])
{
	double sumFit1 = 0, sumFit2 = 0, sumFit3 = 0;

	// sumFit1
	for (int j = 0; j < dimension; j++)
	{
		fv = pow(a[j], 2);
		sumFit1 = sumFit1 + fv;
	}

	// sumFit2
	for (int j = 0; j < dimension; j++)
	{
		fv = 0.5 * (j + 1) * a[j];
		sumFit2 = sumFit2 + fv;
	}
	sumFit2 = pow(sumFit2, 2);

	// sumFit3
	for (int j = 0; j < dimension; j++)
	{
		fv = 0.5 * (j + 1) * a[j];
		sumFit3 = sumFit3 + fv;
	}
	sumFit3 = pow(sumFit3, 4);

	sumFit = sumFit1 + sumFit2 + sumFit3;

	return sumFit;
}

// No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
// Done By: Loh Chia Heung 2301684
double AxisParallel(double a[])
{
	for (int j = 0; j < dimension; j++)
	{
		fv = (j + 1) * pow(a[j], 2);
		sumFit = sumFit + fv;
	}

	return sumFit;
}

// No.6 - Griewank Function +-600
// Done By: Loh Chia Heung 2301684
double Griewank(double a[])
{
	//    double sumFit = 0.0;
	double product = 1.0;

	for (int j = 0; j < dimension; j++)
	{
		sumFit += (a[j] * a[j]) / 4000.0;
		product *= cos(a[j] / sqrt(j + 1));
	}
	return sumFit - product + 1;
}

// No.7 - Sum of Different Powers function +-1.00
// Done By: Yeap Chun Hong 2206352
double SumOfDifferentPowers(double a[])
{
	for (int j = 0; j < dimension; j++)
	{
		fv = pow(fabs(a[j]), j + 2);
		sumFit = sumFit + fv;
	}
	return sumFit;
}

// No.8 - Rotated Hyper-Ellipsoid Function +-65.536
// Done By: Brandon Ting En Junn 2101751
double Rotated(double a[])
{
	sumFit = Sphere(a);
	fv = sumFit;

	for (int i = 0; i < dimension; i++)
	{
		sumFit = sumFit + fv;
	}

	return sumFit;
}

// No.9 - Schwefel 2.22 Function +-5 [Updated]
// Done By: Ling Ji Xiang 2104584
double Schwefel(double a[])
{
	double sum = 0, product = 1.0;
	for (int i = 0; i < dimension; i++)
	{
		double absolute = fabs(a[i]);
		sum += absolute;
		product *= absolute;
	}
	return sum + product;
}

// No.10 - Exponential function Function +-1.00 [f(x) = -1]
// Done By: Ling Ji Xiang 2104584
double Exponential(double a[])
{
	double result = 0.0;

	// calculate the sum of squares
	for (int j = 0; j < dimension; j++)
	{
		fv = pow(a[j], 2);
		sumFit = sumFit + fv;
	}

	// Calculate the exponential of sum of squares
	result = -exp(-0.5 * sumFit);
	return result;
}

/* Benchmark_Range */
void initialiseBenchmark_Range()
{
	switch (BENCHMARK)
	{
	case 1: // No.1 - Sphere Function +-5.12
		benchmarkFunction = "-Sphere";
		rangeMin = rangeMax = 5120;
		break;

	case 2: // No.2 - Ackley Function +-32.768
		benchmarkFunction = "-Ackley";
		rangeMin = rangeMax = 32768;
		break;

	case 3: // No.3 - Rastrigin Function +-5.12
		benchmarkFunction = "-Rastrigin";
		rangeMin = rangeMax = 5120;
		break;

	case 4: // No.4 - Zakharov Function -5, +10
		benchmarkFunction = "-Zakharov";
		rangeMin = 5000;
		rangeMax = 10000;
		break;

	case 5: // No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
		benchmarkFunction = "-AxisParallel";
		rangeMin = rangeMax = 5120;
		break;

	case 6: // No.6 - Griewank Function +-600
		benchmarkFunction = "-Griewank";
		rangeMin = rangeMax = 600000;
		break;

	case 7: // No.7 - Sum of Different Powers function +-1.00
		benchmarkFunction = "-SumOfDifferentPowers";
		rangeMin = rangeMax = 1000;
		break;

	case 8: // No.8 - Rotated Hyper-Ellipsoid Function +-65.536
		benchmarkFunction = "-Rotated";
		rangeMin = rangeMax = 65536;
		break;

	case 9: // No.9 - Schwefel 2.22 Function +-5 [Updated]
		benchmarkFunction = "-Schwefel";
		rangeMin = rangeMax = 5000;
		break;

	case 10: // No.10 - Exponential function Function +-1.00 [f(x) = -1]
		benchmarkFunction = "-Exponential";
		rangeMin = rangeMax = 1000;
		break;

	default: // ERROR
		cout << "Error BENCHMARK: Invalid Benchmark Function\n\nPress Any Key to Exit..." << endl;
		getch();
		exit(0);
	}
}

/* Fitness Function */
double Fitness(double a[])
{
	switch (BENCHMARK)
	{
	case 1: // No.1 - Sphere Function +-5.12
		return Sphere(a);

	case 2: // No.2 - Ackley Function +-32.768
		return Ackley(a);

	case 3: // No.3 - Rastrigin Function +-5.12
		return Rastrigin(a);

	case 4: // No.4 - Zakharov Function -5, +10
		return Zakharov(a);

	case 5: // No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
		return AxisParallel(a);

	case 6: // No.6 - Griewank Function +-600
		return Griewank(a);

	case 7: // No.7 - Sum of Different Powers function +-1.00
		return SumOfDifferentPowers(a);

	case 8: // No.8 - Rotated Hyper-Ellipsoid Function +-65.536
		return Rotated(a);

	case 9: // No.9 - Schwefel 2.22 Function +-5 [Updated]
		return Schwefel(a);

	case 10: // No.10 - Exponential function Function +-1.00 [f(x) = -1]
		return Exponential(a);

	default: // ERROR
		cout << "Error BENCHMARK: Invalid Benchmark Function\n\nPress Any Key to Exit..." << endl;
		getch();
		exit(0);
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// External Function
//------------------------------------------------------------------------------------------------------------------------------
int charToInt(char c)
{
	if (c >= '0' && c <= '9')
	{
		return c - '0';
	}
	else
	{
		// ERROR
		cout << "Error charToInt(): Invalid Conversion\n\nPress Any Key to Exit..." << endl;
		getch();
		exit(0);
	}
}

string addBinary(string a, string b)
{
	string result = ""; 	 // Initialize the result as an empty string
	int s = 0;				 // Initialize the carry

	// Ensure both strings are of the same length by padding with leading zeros
	int n = max(a.size(), b.size());
	while (a.size() < n)
		a.insert(a.begin(), '0');
	while (b.size() < n)
		b.insert(b.begin(), '0');

	// Traverse both strings from right to left
	for (int i = n - 1; i >= 0; i--)
	{
		int sum = (a[i] - '0') + (b[i] - '0') + s;		// Calculate the sum of the current digits and carry
		result.insert(result.begin(), (sum % 2) + '0'); // Insert the current bit to the result
		s = sum / 2;									// Calculate the new carry
	}

	// If there's a carry left, add it to the result
	if (s != 0)
	{
		result.insert(result.begin(), '1');
	}

	return result;
}

void GA_TO_GA_COMBINATION()
{
	int index = 0;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < TECHNIQUE; j++)
		{
			GA_COMBINATION[i][j] = charToInt(GA[index++]);
		}
	}
}

bool isValidCombination(int combination[4][TECHNIQUE])
{
	int zeros = 0, ones = 0;

	// Check Each Operation Technique
	for (int i = 0; i < 4; i++)
	{

		zeros = ones = 0;

		// Check Each Technique
		for (int j = 0; j < TECHNIQUE; j++)
		{
			if (combination[i][j] == 0)
			{
				zeros++;
			}
			if (combination[i][j] == 1)
			{
				ones++;
			}
		}

		if (zeros == 0 || ones == 0)
		{
			return false;
		}
	}

	return true;
}

void resetExperiment()
{
	for (int i = 0; i < pSize; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			chromosome[i][j] = 0;

			if (i < 4)
			{
				paroff[i][j] = 0;
				tfit[i] = 0;
			}

			if (i == 0)
			{
				lowestGene[j] = 0;
			}
		}

		fit[i] = 0;
	}

	r, gcp, gmp = 0;
	crb, mb1, mb2 = 0;
	mb1v, mb2v = 0;
	fv, sumFit = 0;
	fit1, fit2 = 0;

	lFvIndex = 0;
	lFv = 0;

	lowestGeneFV = 0;
	dynamicTournamentSize = tournamentSize;
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Selection Function
// Done By: Yeap Chun Hong 2206352
//------------------------------------------------------------------------------------------------------------------------------
int RouletteWheelSelection2(double fitness[])
{
    
	double inverseFit[pSize];
	double totalFit = 0;

	for (int i = 0; i < pSize; i++)
	{
		inverseFit[i] = 1 / fitness[i]; //inverse so that smaller fitness value will have a higher portion
		totalFit += inverseFit[i];
		//cout <<fitness[i] <<"\t" <<inverseFit[i] <<endl;
	}
	//cout <<totalFit <<endl;

	//calculate the cumulative probability and normalise it to [0,1]
	double cumulativeProbability[pSize];
	cumulativeProbability[0] = inverseFit[0] / totalFit;
	//cout <<"cumulativeProbability" << cumulativeProbability[0] <<endl;
	for (int i = 1; i < pSize; i++)
	{
		cumulativeProbability[i] = cumulativeProbability[i - 1] + (inverseFit[i] / totalFit);
		//cout <<i << "\t"<<cumulativeProbability[i] <<endl;
	}

	// generate number from 0 to 1
	double spin = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);	//cout <<"spin: " <<spin <<endl;
	//cout << spin << endl;

	for (int i = 0; i < pSize; i++)
	{
		if (spin <= cumulativeProbability[i])
		{
			//cout << "chosen: " << i << endl;
			return i;
		}
	}

	// In case of rounding errors, return the last individual
	return pSize;
}

void RouletteWheelSelection(double fitness[], int &parent1, int &parent2)
{
	parent1 = RouletteWheelSelection2(fitness);
	parent2 = RouletteWheelSelection2(fitness);

	// cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;
	// getch();
}

int TournamentSelection2(double fitness[], int tournamentSize)
{
	int best = -1;
	double bestFitness = numeric_limits<double>::max();
	bool selectedIndices[pSize] = { false };

	//Randomly select unique individuals and perform tournament  
	for (int i = 0; i < tournamentSize; i++)
	{
		int index;
		do
		{
			index = rand() % pSize;
		} while (selectedIndices[index]);
		//cout <<"index "<< index <<endl;

		selectedIndices[index] = true;

		//cout <<"fitness[index] "<< fitness[index]<<" bestFitness "<< bestFitness <<endl;
		if (fitness[index] < bestFitness)
		{
			best = index;
			bestFitness = fitness[index];
			//cout <<"best "<< best<<" bestFitness "<< bestFitness <<endl;
		}
	}
	//getch();
	return best;
}

void TournamentSelection(double fitness[], int tournamentSize, int &parent1, int &parent2)
{
	// Select first parent
	parent1 = TournamentSelection2(fitness, tournamentSize);

	// Select second parent
	parent2 = TournamentSelection2(fitness, tournamentSize);

	// cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;
	// getch();
}

int LinearRankingSelection2(double fitness[])
{
	//selection pressure
	double max = 1.1;
	double min = 2 - max;


	//sort fitness value in descending order
	sort(fitness, fitness + pSize, greater<double>());
	double probability[pSize];
	double totalFit = 0;
	for (int i = 0; i < pSize; i++)
	{
		probability[i] = (min + (max - min) * i / (pSize - 1)) / pSize;
		totalFit += probability[i];

	}

	//calculate the cumulative probability and normalise it to [0,1]
	double cumulativeProbability[pSize];
	cumulativeProbability[0] = probability[0] / totalFit;
	//cout <<"i" <<"\t"<< "Fitness"<<"\t"<< "Probability"<< "\t" << "Cumulative Prob" << endl;
	//cout <<"0" <<"\t"<< fitness[0]<<"\t"<< probability[0]<< "\t" <<cumulativeProbability[0] << endl;

	//cout <<"cumulativeProbability" << cumulativeProbability[0] <<endl;
	for (int i = 1; i < pSize; i++)
	{
		cumulativeProbability[i] = cumulativeProbability[i - 1] + (probability[i] / totalFit);
		//cout <<i <<"\t"<< fitness[i]<<"\t"<< probability[i]<< "\t" <<cumulativeProbability[i] << endl;
		//cout <<i << "\t"<<cumulativeProbability[i] <<endl;
	}

	// generate number from 0 to 1
	double spin = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);	//cout <<"spin: " <<spin <<endl;
	//cout << spin <<endl;

	for (int i = 0; i < pSize; i++)
	{
		if (spin <= cumulativeProbability[i])
		{
			//cout << "chosen: "<<i<<endl;
			//getch();
			return i;
		}
	}

	// In case of rounding errors, return the last individual
	return pSize;
}

void LinearRankingSelection(double fitness[], int& parent1, int& parent2)
{
	// Select first parent
	parent1 = LinearRankingSelection2(fitness);

	// Select second parent
	parent2 = LinearRankingSelection2(fitness);

	//getch();
}

double CalculateDiversity(double fitness[])
{
	
	double mean = 0.0;
	double sumSquare = 0.0;
	double sumFitness =0.0;
	for (int i =0; i < pSize; i++){
		sumFitness += fitness[i];
	}
	//calculate mean of fitness
	for (int i =0; i < pSize; i++){
		mean += fitness[i]/sumFitness;
	}
	
	mean = mean/pSize;
	
	//calculate sum of square
	for (int i=0; i < pSize; i++){
		sumSquare += (fitness[i]/sumFitness - mean) * (fitness[i]/sumFitness - mean);
	}
	
	double stdDev = sqrt(sumSquare/pSize);
	return stdDev;
}

int AdjustTournamentSize(double diversity, int currentSize, int minSize, int maxSize) {
	//cout <<"Diversity: "<< diversity <<"\t" <<"Size: " <<currentSize<<endl;
	//getch();
	
    if (diversity > highDiverseThreshold) {  // if diversity > 0.04, lower tournamentsize
        return max(minSize, currentSize - 1);
    } else if (diversity < lowDiverseThreshold) {  // if diversity < 0.01, lower tournamentsize
        return min(maxSize, currentSize + 1);
    }
    
    //remain the same if 0.01 < diversity < 0.04
    return currentSize;
}

void DynamicTournamentSelection(double fitness[], int& dynamicTournamentSize, int& parent1, int& parent2)
{
	double diversity = CalculateDiversity(fitness);
	dynamicTournamentSize = AdjustTournamentSize(diversity,dynamicTournamentSize,3,7);	
	// Select first parent
    parent1 = TournamentSelection2(fitness, dynamicTournamentSize);

    // Select second parent
    parent2 = TournamentSelection2(fitness, dynamicTournamentSize);
    
    //getch();
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Crossover Function
// Done By: Loh Chia Heung 2301684
//------------------------------------------------------------------------------------------------------------------------------
void UniformCrossover(double chromosome[][dimension], double paroff[][dimension], double dcp, int p1, int p2)
{
	double gcp = (rand() % 1000);
	gcp = gcp / 1000;

	if (gcp <= dcp)
	{
		// Perform Uniform Crossover
		for (int j = 0; j < dimension; j++)
		{
			if (rand() % 2 == 0)
			{
				// 0 --> Gene is copied from Parent 1
				paroff[2][j] = chromosome[p1][j];
				paroff[3][j] = chromosome[p2][j];
			}
			else
			{
				// 1 --> Gene is copied from Parent 2
				paroff[2][j] = chromosome[p2][j];
				paroff[3][j] = chromosome[p1][j];
			}
		}
	}
	else
	{
		// No crossover, directly copy parents to offspring
		for (int j = 0; j < dimension; j++)
		{
			paroff[2][j] = chromosome[p1][j]; // Offspring1 --> Parent 1
			paroff[3][j] = chromosome[p2][j]; // Offspring2 --> Parent 2
		}
	}
}

void ShuffleCrossover(double chromosome[][dimension], double paroff[][dimension], double dcp, int p1, int p2)
{

	double gcp = (rand() % 1000);
	gcp = gcp / 1000;

	if (gcp <= dcp)
	{
		double temp[dimension];

		// Shuffle chromosomes between parents
		for (int j = 0; j < dimension; j++)
		{
			temp[j] = (rand() % 2 == 0) ? chromosome[p1][j] : chromosome[p2][j];
		}

		// Assign shuffled chromosomes to offspring
		for (int j = 0; j < dimension; j++)
		{
			paroff[2][j] = temp[j];
			paroff[3][j] = (rand() % 2 == 0) ? chromosome[p1][j] : chromosome[p2][j];
		}
	}
	else
	{
		// No crossover --> offspring are exact copies of parents
		for (int j = 0; j < dimension; j++)
		{
			paroff[2][j] = chromosome[p1][j]; // Offspring1 --> Parent 1
			paroff[3][j] = chromosome[p2][j]; // Offspring2 --> Parent 2
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Mutation Function
// Done By: Brandon Ting En Junn 2101751
//------------------------------------------------------------------------------------------------------------------------------
void ReversingMutation()
{
	gmp = (rand() % 1000000);
	gmp = gmp / 1000000;

	// Parent 1
	if (gmp <= dmp)
	{
		mb1 = getrandom(0, dimension - 1);

		mb1v = paroff[2][mb1];

		if (mb1 == 0)
		{
			// Parent 1
			paroff[2][mb1] = paroff[2][dimension - 1];
			paroff[2][dimension - 1] = mb1v;
		}
		else
		{
			// Parent 1
			paroff[2][mb1] = paroff[2][mb1 - 1];
			paroff[2][mb1 - 1] = mb1v;
		}
	}

	gmp = (rand() % 1000000);
	gmp = gmp / 1000000;

	// Parent 2
	if (gmp <= dmp)
	{
		mb2 = getrandom(0, dimension - 1);

		mb2v = paroff[2][mb2];

		if (mb2 == 0)
		{
			// Parent 1
			paroff[3][mb2] = paroff[3][dimension - 1];
			paroff[3][dimension - 1] = mb2v;
		}
		else
		{
			// Parent 1
			paroff[3][mb1] = paroff[3][mb2 - 1];
			paroff[3][mb1 - 1] = mb2v;
		}
	}
}

void RandomMutation()
{
	// Parent 1 & Parent 2
	for (int i = 2; i < 4; i++)
	{
		gmp = (rand() % 1000000);
		gmp = gmp / 1000000;

		mb1 = getrandom(0, dimension - 1);
		mb2 = getrandom(0, dimension - 1);

		if (gmp <= dmp)
		{
			r = getrandom(-rangeMin, rangeMax);
			r = r / rangeDiv;
			paroff[i][mb1] = r;
		}

		gmp = (rand() % 1000000);
		gmp = gmp / 1000000;

		if (gmp <= dmp)
		{
			r = getrandom(-rangeMin, rangeMax);
			r = r / rangeDiv;
			paroff[i][mb2] = r;
		}
	}
}

void SimpleInversionMutation()
{
	for (int i = 2; i < 4; i++)
	{
		gmp = (rand() % 1000000);
		gmp = gmp / 1000000;

		if (gmp <= dmp)
		{
			do {
				mb1 = getrandom(0, dimension - 1);
				mb2 = getrandom(0, dimension - 1);

				if (mb1 > mb2) {
					double temp = mb1;
					mb1 = mb2;
					mb2 = temp;
				}
			} while (mb1 == mb2);

			int iterations = (mb2 - mb1 + 1) / 2;
			for (int j = 0; j < iterations; j++, mb1++, mb2--)
			{
				double temp = paroff[i][mb1];
				paroff[i][mb1] = paroff[i][mb2];
				paroff[i][mb2] = temp;
			}
		}
	}
}

void ImitatingMutation()
{
	/*
		The Imitating Mutation selects a random gene to be mutated.
		Then, the genes are iterated to find the smallest/largest/optimum gene in the chromosome to be imitated.
		The selected random gene is replaced with the optimum gene.
	*/
	// Child 1 & Child 2 Loop
	for (int i = 2; i < 4; i++)
	{
		gmp = (rand() % 1000000);
		gmp = gmp / 1000000;

		// Mutation Occurs
		if (gmp <= dmp)
		{
			double optimumGene = numeric_limits<double>::max();
			mb1 = getrandom(0, dimension - 1);

			// Find optimumGene
			for (int j = 0; j < dimension; j++)
			{
				if (j != mb1) {
					if (paroff[i][j] < optimumGene) {
						optimumGene = paroff[i][j];
					}
				}
			}

			// Gene Mutation
			paroff[i][mb1] = optimumGene;
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Replacement Function
// Done By: Ling Ji Xiang 2104584
//------------------------------------------------------------------------------------------------------------------------------
void WeakParentReplacement(double chromosome[pSize][dimension], double fit[pSize], double paroff[4][dimension], double tfit[4], int parent1, int parent2)
{
	int chosenParent, notChosenParent;
	int chosenChild, notChosenChild;
	
	if (fit[parent1] > fit[parent2]) //choose weak parent
	{
		chosenParent = parent1; 
		notChosenParent = parent2; 
	}
	else
	{
		chosenParent = parent2; 
		notChosenParent = parent1;
	}
	
	if (tfit[2] < tfit[3]) //choose strong child
	{
		chosenChild = 2; 
		notChosenChild = 3; 
	}
	else
	{
		chosenChild = 3; 
		notChosenChild = 2; 
	}
	
	if (fit[chosenParent] > tfit[chosenChild])
	{
		for (int i = 0; i < dimension; i++)
		{
			chromosome[chosenParent][i] = paroff[chosenChild][i];
		}
		fit[chosenParent] = tfit[chosenChild];
	}
	
	if (fit[notChosenParent] > tfit[notChosenChild])
	{
		for (int i = 0; i < dimension; i++)
		{
			chromosome[notChosenParent][i] = paroff[notChosenChild][i];
		}
		fit[notChosenParent] = tfit[notChosenChild];
	}
}

void BothParentReplacement(double chromosome[pSize][dimension], double fit[pSize], double paroff[4][dimension], double tfit[4], int parent1, int parent2)
{
	for (int i = 0; i < dimension; i++)
	{
		chromosome[parent1][i] = paroff[2][i]; // Replace genes of parent1 with offspring 1
		chromosome[parent2][i] = paroff[3][i]; // Replace genes of parent2 with offspring 2
	}
	fit[parent1] = tfit[2]; // Update fitness of parent1
	fit[parent2] = tfit[3]; // Update fitness of parent2
}
//------------------------------------------------------------------------------------------------------------------------------

int main()
{
#ifdef EXIT // ERROR
	cout << "Error MINI_PROJECT: Invalid Parameter Settings\n\nPress Any Key to Exit..." << endl;
	getch();
	exit(0);
#endif

	srand(time(0));

	//---------------------------------------------------------------------------------------------------------------------------
	// GA Automation
	// Done By: Brandon Ting En Junn 2101751
	//---------------------------------------------------------------------------------------------------------------------------
	int GACounter = 1, GAModel = 0;
	const int GALength = pow(2, GA.length());

	do
	{
		//Check GA_COMBINATION e.g. (S)00 (C)00 (M)00 (R)00
		GA_TO_GA_COMBINATION();
		if (!isValidCombination(GA_COMBINATION))
		{
			GA = addBinary(GA, "1");
			continue; //[Skip] e.g. (S)00 (C)00 (M)00 (R)00
		}

		//[Proceed] e.g. (S)01 (C)01 (M)01 (R)01
		GAModel++;
		//GA Directory e.g. GA1
		string gaDirectory = "GA" + to_string(GAModel);
		CreateDirectory(gaDirectory.c_str(), NULL);									
		string gaInfo = gaDirectory + "\\\\" + gaDirectory + ".txt";	//GA Info e.g. GA1.txt
		ofstream outfileo1Info(gaInfo.c_str(), ios::trunc);

		cout << gaDirectory << endl;
		outfileo1Info << gaDirectory << endl << endl;

		//---------------------------------------------------------------------------------------------------------------------------
		// Benchmark Automation
		// Done By: Brandon Ting En Junn 2101751
		//---------------------------------------------------------------------------------------------------------------------------
		for (; BENCHMARK < 11; BENCHMARK++)
		{
			initialiseBenchmark_Range(); //Set benchmark mode
			//Benchmark Directory e.g. GA1/GA1-Sphere
			string benchmarkDirectory = gaDirectory + "\\\\" + gaDirectory + benchmarkFunction;
			CreateDirectory(benchmarkDirectory.c_str(), NULL);

			//---------------------------------------------------------------------------------------------------------------------------
			// Experiment Automation
			// Done By: Brandon Ting En Junn 2101751
			//---------------------------------------------------------------------------------------------------------------------------
			for (int experiment = 0; experiment < 10; experiment++)
			{
				// Result File e.g. GA1/GA1-Sphere/GA1-Sphere-Result1.txt
				string outfile1 = benchmarkDirectory + "\\\\" + gaDirectory + benchmarkFunction + "-Result" + to_string(experiment + 1) + ".txt";
				ofstream outfileo1(outfile1.c_str(), ios::trunc);

				cout << "Experiment " << experiment + 1 << "..." << endl;

				// CPU Time
				clock_t start, end;
				start = clock();

				//---------------------------------------------------------------------------------------------------------------------------
				// Generate Population
				//---------------------------------------------------------------------------------------------------------------------------
#if MINI_PROJECT == -1
				cout << "Generate Population..." << endl;
#endif

				for (int i = 0; i < pSize; i++)
				{
					for (int j = 0; j < dimension; j++)
					{
						r = getrandom(-rangeMin, rangeMax);
						r = r / rangeDiv;
						chromosome[i][j] = r;
					}

#if MINI_PROJECT == -1
					if (i == 0 || i == pSize - 1)
					{
						cout << "Chromosome " << i + 1 << endl;
						for (int j = 0; j < dimension; j++)
						{
							cout << setprecision(6) << chromosome[i][j] << "\t";
						}
						cout << endl
							 << endl;
					}
#endif

				}
				//---------------------------------------------------------------------------------------------------------------------------

				//---------------------------------------------------------------------------------------------------------------------------
				// Fitness Evaluation 1
				//---------------------------------------------------------------------------------------------------------------------------
				for (int i = 0; i < pSize; i++)
				{
					fit[i] = Fitness(chromosome[i]);
					//cout<<setprecision(6)<<fit[i]<<endl;
					sumFit = 0;
				}

				lFv = numeric_limits<double>::max();

				for (int i = 0; i < pSize; i++)
				{
					if (fit[i] < lFv)
					{
						lFv = fit[i];
						lFvIndex = i;
					}
				}

				lowestGeneFV = lFv;
				//cout<<"The best is chromosome "<<lFvIndex+1<<" with fitness of "<<lowestGeneFV<<endl;
				for (int j = 0; j < dimension; j++)
				{
					lowestGene[j] = chromosome[lFvIndex][j];
					//cout<<lowestGene[j]<<" ";
				}
				//cout<<endl;
				//outfileo1<<"Gen\tMinimum"<<endl;
				//---------------------------------------------------------------------------------------------------------------------------

				//---------------------------------------------------------------------------------------------------------------------------
				// Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
				//---------------------------------------------------------------------------------------------------------------------------
				for (int i = 0; i < gen; i++)
				{
					//------------------------------------------------------------------------------------------------------------------------
					// Selection Operation
					// Done By: Yeap Chun Hong 2206352
					//------------------------------------------------------------------------------------------------------------------------
#if MINI_PROJECT == -1
					cout << "Selection Operation..." << endl;
#endif

					int parent1 = 0, parent2 = 0;

					//************************************************************************************************************************
#if MINI_PROJECT == -1
					/* Best Selection */
					if (GA_COMBINATION[0][1] == 1)
					{
						RouletteWheelSelection(fit, parent1, parent2);
						cout << "Roulette Wheel Selection" << endl;
						if (i == 0) { outfileo1Info << "S: Roulette Wheel Selection" << endl; }
					}
#elif MINI_PROJECT == -2
					/* Manual Selection */
					if (GA_COMBINATION[0][1] == 1)
					{
						RouletteWheelSelection(fit, parent1, parent2);
						if (i == 0) { outfileo1Info << "S: Roulette Wheel Selection" << endl; }
					}
#else
					/* Roulette Wheel Selection */
					if (GA_COMBINATION[0][0] == 1)
					{
						RouletteWheelSelection(fit, parent1, parent2);
						if (i == 0) { outfileo1Info << "S: Roulette Wheel Selection" << endl; }
					}

					/* Tournament Selection */
					if (GA_COMBINATION[0][1] == 1)
					{
						TournamentSelection(fit, tournamentSize, parent1, parent2);
						if (i == 0) { outfileo1Info << "S: Tournament Selection" << endl; }
					}
#if MINI_PROJECT == 1
					/* Linear Ranking Selection */
					if (GA_COMBINATION[0][2] == 1)
					{
						LinearRankingSelection(fit, parent1, parent2);
						if (i == 0) { outfileo1Info << "S: Linear Ranking Selection" << endl; }
					}

					/* Dynamic Tournament Selection */
					if (GA_COMBINATION[0][3] == 1)
					{
						DynamicTournamentSelection(fit, dynamicTournamentSize, parent1, parent2);
						if (i == 0) { outfileo1Info << "S: Dynamic Tournament Selection" << endl; }
					}
#endif
#endif
					//************************************************************************************************************************

					tfit[0] = fit[parent1];
					tfit[1] = fit[parent2];

					for (int j = 0; j < dimension; j++)
					{
						paroff[0][j] = chromosome[parent1][j];
						paroff[1][j] = chromosome[parent2][j];
					}

#if MINI_PROJECT == -1
					// Selected Parent 1
					cout << "Chromosome " << parent1 + 1 << endl;
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << chromosome[parent1][j] << "\t";
					}
					cout << endl << endl;

					// Selected Parent 2
					cout << "Chromosome " << parent2 + 1 << endl;if (i == 0) {  }
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << chromosome[parent2][j] << "\t";
					}
					cout << endl << endl;
#endif
					//------------------------------------------------------------------------------------------------------------------------

					//------------------------------------------------------------------------------------------------------------------------
					// Crossover Operation
					// Done By: Loh Chia Heung 2301684
					//------------------------------------------------------------------------------------------------------------------------
#if MINI_PROJECT == -1
					cout << "Crossover Operation..." << endl;
#endif

					//************************************************************************************************************************
#if MINI_PROJECT == -1
					/* Best Crossover */
					if (GA_COMBINATION[1][1] == 1)
					{
						UniformCrossover(chromosome, paroff, dcp, parent1, parent2);
						cout << "Uniform Crossover" << endl;
						if (i == 0) { outfileo1Info << "C: Uniform Crossover" << endl; }
					}
#elif MINI_PROJECT == -2
					/* Manual Crossover */
					if (GA_COMBINATION[1][1] == 1)
					{
						UniformCrossover(chromosome, paroff, dcp, parent1, parent2);
						if (i == 0) { outfileo1Info << "C: Uniform Crossover" << endl; }
					}
#else
					/* Uniform Crossover */
					if (GA_COMBINATION[1][0] == 1)
					{
						UniformCrossover(chromosome, paroff, dcp, parent1, parent2);
						if (i == 0) { outfileo1Info << "C: Uniform Crossover" << endl; }
					}

					/* Shuffle Crossover */
					if (GA_COMBINATION[1][1] == 1)
					{
						ShuffleCrossover(chromosome, paroff, dcp, parent1, parent2);
						if (i == 0) { outfileo1Info << "C: Shuffle Crossover" << endl; }
					}
#if MINI_PROJECT == 1
					/* XXX Crossover */
					//if (GA_COMBINATION[1][2] == 1)
					//{
					//	
					//	if (i == 0) { outfileo1Info << "C: XXX Crossover" << endl; }
					//}

					/* XXX Crossover */
					//if (GA_COMBINATION[1][3] == 1)
					//{
					//	
					//	if (i == 0) { outfileo1Info << "C: XXX Crossover" << endl; }
					//}
#endif
#endif
					//************************************************************************************************************************

#if MINI_PROJECT == -1
					// Crossover Child 1
					cout << "Child 1" << endl;
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << paroff[2][j] << "\t";
					}
					cout << endl << endl;

					// Crossover Child 2
					cout << "Child 2" << endl;
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << paroff[3][j] << "\t";
					}
					cout << endl << endl;
#endif
					//------------------------------------------------------------------------------------------------------------------------

					//------------------------------------------------------------------------------------------------------------------------
					// Mutation Operation
					// Done By: Brandon Ting En Junn 2101751
					//------------------------------------------------------------------------------------------------------------------------
#if MINI_PROJECT == -1
					cout << "Mutation Operation..." << endl;
#endif

					//************************************************************************************************************************
#if MINI_PROJECT == -1
					/* Best Mutation */
					if (GA_COMBINATION[2][1] == 1)
					{
						ReversingMutation();
						cout << "Reversing Mutation" << endl;
						if (i == 0) { outfileo1Info << "M: Reversing Mutation" << endl; }
					}
#elif MINI_PROJECT == -2
					/* Manual Mutation */
					if (GA_COMBINATION[2][1] == 1)
					{
						ImitatingMutation();
						if (i == 0) { outfileo1Info << "M: Imitating Mutation" << endl; }
					}
#else
					/* Reversing Mutation */
					if (GA_COMBINATION[2][0] == 1)
					{
						ReversingMutation();
						if (i == 0) { outfileo1Info << "M: Reversing Mutation" << endl; }
					}

					/* Random Mutation */
					if (GA_COMBINATION[2][1] == 1)
					{
						RandomMutation();
						if (i == 0) { outfileo1Info << "M: Random Mutation" << endl; }
					}
#if MINI_PROJECT == 1
					/* Simple Inversion Mutation */
					if (GA_COMBINATION[2][2] == 1)
					{
						SimpleInversionMutation();
						if (i == 0) { outfileo1Info << "M: Simple Inversion Mutation" << endl; }
					}

					/* Imitating Mutation */
					if (GA_COMBINATION[2][3] == 1)
					{
						ImitatingMutation();
						if (i == 0) { outfileo1Info << "M: Imitating Mutation" << endl; }
					}
#endif
#endif
					//************************************************************************************************************************

#if MINI_PROJECT == -1
					// Mutation Child 1
					cout << "Child 1" << endl;
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << paroff[2][j] << "\t";
					}
					cout << endl << endl;

					// Mutation Child 2
					cout << "Child 2" << endl;
					for (int j = 0; j < dimension; j++)
					{
						cout << setprecision(6) << paroff[3][j] << "\t";
					}
					cout << endl << endl;
#endif
					//------------------------------------------------------------------------------------------------------------------------

					//------------------------------------------------------------------------------------------------------------------------
					// Fitness Evaluation 2
					//------------------------------------------------------------------------------------------------------------------------
					fit1 = Fitness(paroff[2]);
					sumFit = 0;

					fit2 = Fitness(paroff[3]);
					sumFit = 0;

					tfit[2] = fit1;
					tfit[3] = fit2;
					//------------------------------------------------------------------------------------------------------------------------

					//------------------------------------------------------------------------------------------------------------------------
					// Replacement Operation
					// Done By: Ling Ji Xiang 2104584
					//------------------------------------------------------------------------------------------------------------------------
#if MINI_PROJECT == -1
					cout << "Mutation Operation..." << endl;
#endif

					//************************************************************************************************************************
#if MINI_PROJECT == -1
					/* Best Replacement */
					if (GA_COMBINATION[3][1] == 1)
					{
						WeakParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
						cout << "Weak Parent Replacement" << endl << endl << endl;
						if (i == 0) { outfileo1Info << "R: Weak Parent Replacement" << endl << endl << endl; }
					}
#elif MINI_PROJECT == -2
					/* Manual Replacement */
					if (GA_COMBINATION[3][1] == 1)
					{
						WeakParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
						if (i == 0) { outfileo1Info << "R: Weak Parent Replacement" << endl << endl << endl; }
					}
#else
					/* Weak Parent Replacement */
					if (GA_COMBINATION[3][0] == 1)
					{
						WeakParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
						if (i == 0) { outfileo1Info << "R: Weak Parent Replacement" << endl << endl << endl; }
					}

					/* Both Parent Replacement */
					if (GA_COMBINATION[3][1] == 1)
					{
						BothParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
						if (i == 0) { outfileo1Info << "R: Both Parent Replacement" << endl << endl << endl; }
					}
#if MINI_PROJECT == 1
					/* XXX Replacement */
					//if (GA_COMBINATION[3][2] == 1)
					//{
					//	
					//	if (i == 0) { outfileo1Info << "R: XXX Replacement" << endl; }
					//}

					/* XXX Replacement */
					//if (GA_COMBINATION[3][3] == 1)
					//{
					//	
					//	if (i == 0) { outfileo1Info << "R: XXX Replacement" << endl; }
					//}
#endif
#endif
					//************************************************************************************************************************
					//------------------------------------------------------------------------------------------------------------------------

#if MINI_PROJECT == -1
					cout << "SUCCESSFULLY COMPLETED DEMO...\nPress Any Key to Exit..." << endl;
					getch();
					exit(0);
#endif

					lFv = numeric_limits<double>::max();
					for (int j = 0; j < pSize; j++)
					{
						if (fit[j] < lFv)
						{
							lFv = fit[j];
							lFvIndex = j;
						}

						//cout<<lFv<<" "<<lowestGeneFV<<endl;
						if (lFv < lowestGeneFV)
						{
							lowestGeneFV = lFv;

							for (int k = 0; k < dimension; k++)
							{
								lowestGene[k] = chromosome[lFvIndex][k];
								//cout<<lowestGene[k]<<" ";
							}
						}
					}

					fit1 = 0;
					fit2 = 0;

					outfileo1 << setprecision(6) << lFv << endl;
				} // Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)] LOOP
				//------------------------------------------------------------------------------------------------------------------------

				lFv = numeric_limits<double>::max();
				for (int j = 0; j < pSize; j++)
				{
					if (fit[j] < lFv)
					{
						lFv = fit[j];
						lFvIndex = j;
					}
				}

				outfileo1 << endl << endl;
				cout << "Result" << endl;

				// Output Final Generation Lowest
				//cout<<setprecision(6)<<lFv<<" "<<lFvIndex+1<<endl<<endl;
				//for(int j = 0 ; j < dimension ; j++)
				//{
				//	cout<<setprecision(6)<<chromosome[lFvIndex][j]<<"\t";
				//	outfileo1<<setprecision(6)<<chromosome[lFvIndex][j]<<"\n";
				//}

				// Output Lowest Gene in Experiment
				cout << setprecision(6) << lowestGeneFV << endl << endl;
				for (int j = 0; j < dimension; j++)
				{
					cout << setprecision(6) << lowestGene[j] << "\t";
					outfileo1 << setprecision(6) << lowestGene[j] << "\n";
				}

				cout << endl;
				outfileo1 << endl;

				end = clock();
				resetExperiment();
				cout << "Time required for execution: " << (double)(end - start) / CLOCKS_PER_SEC << " seconds." << "\n\n";
				cout << endl << endl;
				outfileo1 << (double)(end - start) / CLOCKS_PER_SEC << "\n\n";
				outfileo1.close();
			} //Experiment Automation LOOP
			//------------------------------------------------------------------------------------------------------------------------

		} //Benchmark Automation LOOP
		//------------------------------------------------------------------------------------------------------------------------

#if MINI_PROJECT == -2
		cout << "SUCCESSFULLY COMPLETED...\nPress Any Key to Exit..." << endl;
		getch();
		exit(0);
#endif

		outfileo1Info.close();
		BENCHMARK = 1;			 //Reset benchmark mode
		GA = addBinary(GA, "1"); //Next GA combination
	} while (++GACounter <= GALength); //GA Automation LOOP
	//------------------------------------------------------------------------------------------------------------------------

	cout << "SUCCESSFULLY COMPLETED...\nPress Any Key to Exit..." << endl;
	getch();

	return 0;
}
