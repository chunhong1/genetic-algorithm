#define _USE_MATH_DEFINES
#include <iostream>
#include <conio.h>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <windows.h>
using namespace std;

/*************************************** CODE GUIDELINE ***************************************/
//Parameter Settings

//Benchmark Function
//Selection Function
//Crossover Function
//Mutation Function
//Replacement Function

//External Function

//Experiment
	//Generate Population
	//Fitness Evaluation 1
	//Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
		//Selection Operation
		//Crossover Operation
		//Mutation Operation
		//Fitness Evaluation 2
		//Replacement Operation
/*************************************** CODE GUIDELINE ***************************************/

//------------------------------------------------------------------------------------------------------------------------------
//Parameter Settings
//------------------------------------------------------------------------------------------------------------------------------
#define DEMO 0							//Generation = 1 + CMD Output
#define EXPERIMENT 10					//No. of Experiments
const int BENCHMARK = 1;				//Benchmark Function (1 - 10)

//Selection (Roulette Wheel)
//Roulette Wheel
//Tournament

//Crossover (Uniform)
//Uniform
//Shuffle

//Mutation (Reversing)
//Reversing
//Random

//Replacement (Weak Parent)
//Weak Parent
//Both Parent
//------------------------------------------------------------------------------------------------------------------------------

#define getrandom(min,max) (static_cast<long long>(rand()) * (max - min + 1) / RAND_MAX) + min
#if DEMO == 1
	#define gen 1
#else
	#define gen 2000        			//number of iterations (number of generations)
#endif
#define pSize 40           				//number of chromosomes (population size)
#define dimension 30       				//number of bits (dimension size)

float chromosome[pSize][dimension]; 	//chromosome
float paroff[4][dimension];         	//parent and offspring
float fit[pSize];                   	//fitness value for each chromosome
float r = 0, dcp = 0.7, dmp = 0.01, gcp = 0, gmp = 0;
int crb = 0 , mb1 = 0 , mb2 = 0;
int rp1 = 0, rp2 = 0;
float mb1v = 0 , mb2v = 0;
float fv = 0, sumFit=0;
float fit1 = 0 , fit2 = 0;
float tfit[4];

int lFvIndex = 0;
float lFv = 0;

double lowestGene[dimension];
double lowestGeneFV = 0;

const float pi = M_PI;

//------------------------------------------------------------------------------------------------------------------------------
//Benchmark Function
//------------------------------------------------------------------------------------------------------------------------------
//No.1 - Sphere Function +-5.12
float Sphere(float a[])
{
	for(int j = 0 ; j < dimension ; j++) 
	{
		fv = pow(a[j],2);
		sumFit = sumFit + fv;
	}
   
	return sumFit;
}

//No.2 - Ackley Function +-32.768
//Done By: Yeap Chun Hong 2206352
float Ackley (float a[])
{
	float sum1 = 0, sum2 = 0;
	
	for(int j = 0; j < dimension; j++)
	{
		sum1 += pow(a[j], 2);
		sum2 += cos(2 * pi * a[j]);
	}
	
	float term1 = -20 * exp(-0.2 * sqrt(sum1 / dimension));
	float term2 = exp(sum2 / dimension);
	
	sumFit = term1 - term2 + 20 + exp(1);
	return sumFit;
}

//No.3 - Rastrigin Function +-5.12
//Done By: Yeap Chun Hong 2206352
float Rastrigin (float a[])
{
   for(int j = 0 ; j < dimension ; j++) 
   {
	fv = ( pow(a[j],2) ) - ( 10 * cos (2 * pi * a[j]) );
    sumFit = sumFit + fv;
   }
   
   sumFit += 10 * dimension;
   return sumFit;
}

//No.4 - Zakharov Function -5, +10
//Done By: Brandon Ting En Junn 2101751
float Zakharov(float a[])
{
   float sumFit1 = 0, sumFit2 = 0, sumFit3 = 0;
   
   //sumFit1
   for(int j = 0 ; j < dimension ; j++) 
   {
      fv = pow(a[j],2);
      sumFit1 = sumFit1 + fv;
   }
   
   //sumFit2
   for(int j = 0 ; j < dimension ; j++) 
   {
      fv = 0.5 * (j + 1) * a[j];
      sumFit2 = sumFit2 + fv;
   }
   sumFit2 = pow(sumFit2,2);
   
   //sumFit3
   for(int j = 0 ; j < dimension ; j++) 
   {
      fv = 0.5 * (j + 1) * a[j];
      sumFit3 = sumFit3 + fv;
   }
   sumFit3 = pow(sumFit3,4);
   
   sumFit = sumFit1 + sumFit2 + sumFit3;
   
   return sumFit;
}

//No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
//Done By: Loh Chia Heung 2301684
float AxisParallel(float a[])
{  
  
   for(int j = 0; j < dimension; j++) 
   {
      fv = (j + 1) * pow(a[j], 2); 
      sumFit = sumFit + fv;
   }
   
   return sumFit;
}

//No.6 - Griewank Function +-600
//Done By: Loh Chia Heung 2301684
float Griewank(float a[]) {
//    float sumFit = 0.0;
    float product = 1.0;
    
    for (int j = 0; j < dimension; j++) {
        sumFit += (a[j] * a[j]) / 4000.0;
        product *= cos(a[j] / sqrt(j + 1));
    }
    return sumFit - product + 1;
}

//No.7 - Sum of Different Powers function +-1.00
//Done By: Yeap Chun Hong 2206352
float SumOfDifferentPowers(float a[])
{
   for(int j = 0 ; j < dimension ; j++) 
   {
   	fv = pow(fabs(a[j]),j+2);
      sumFit = sumFit + fv;
   }
   return sumFit;
}

//No.8 - Rotated Hyper-Ellipsoid Function +-65.536
//Done By: Brandon Ting En Junn 2101751
float Rotated(float a[])
{
	sumFit = Sphere(a);
	fv = sumFit;
	
	for (int i = 0; i < dimension; i++)
	{
		sumFit = sumFit + fv;
	}
	
	return sumFit;
}

//No.9 - Schwefel 2.22 Function +-5
//Done By: Ling Ji Xiang 2104584
float Schwefel(float a[])
{
	float sum = 0, product = 1.0;
	for(int i = 0; i < dimension; i++)
	{
		float absolute = fabs(a[i]);
		sum += absolute;
		product *= absolute;
	}
	return sum + product;
}

//No.10 - Exponential function Function +-1.00 [f(x) = -1]
//Done By: Ling Ji Xiang 2104584
float Exponential(float a[]) 
{    
    float result = 0.0; 

    //calculate the sum of squares 
    for (int j = 0; j < dimension; j++) 
    { 
        fv = pow(a[j], 2); 
        sumFit = sumFit + fv; 
    } 
    
    //Calculate the exponential of sum of squares 
    result = -exp(-0.5 * sumFit); 
    return result; 
}

/* Benchmark_Range */
string benchmarkFunction = "";
int rangeMin = 0, rangeMax = 0, rangeDiv = 1000;
void initialiseBenchmark_Range()
{
	switch(BENCHMARK)
	{
		case 1: //No.1 - Sphere Function +-5.12
			benchmarkFunction = "-Sphere";
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 2:	//No.2 - Ackley Function +-32.768
			benchmarkFunction = "-Ackley";
		   rangeMin = rangeMax = 32768;
		   break;
		   
		case 3:	//No.3 - Rastrigin Function +-5.12
			benchmarkFunction = "-Rastrigin";
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 4:	//No.4 - Zakharov Function -5, +10
			benchmarkFunction = "-Zakharov";
		   rangeMin = 5000; rangeMax = 10000;
		   break;
		   
		case 5:	//No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
			benchmarkFunction = "-AxisParallel";
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 6:	//No.6 - Griewank Function +-600
			benchmarkFunction = "-Griewank";
		   rangeMin = rangeMax = 600000;
		   break;
		   
		case 7:	//No.7 - Sum of Different Powers function +-1.00
			benchmarkFunction = "-SumOfDifferentPowers";
		   rangeMin = rangeMax = 1000;
		   break;
		   
		case 8: 	//No.8 - Rotated Hyper-Ellipsoid Function +-65.536
			benchmarkFunction = "-Rotated";
			rangeMin = rangeMax = 65536;
		   break;
		   
		case 9: 	//No.9 - Schwefel 2.22 Function +-5 [Updated]
			benchmarkFunction = "-Schwefel";
		   rangeMin = rangeMax = 5000;
		   break;
		   
		case 10: //No.10 - Exponential function Function +-1.00 [f(x) = -1]
			benchmarkFunction = "-Exponential";
		   rangeMin = rangeMax = 1000;
		   break;
		   
		default: //ERROR
		   cout << "Error: Invalid Benchmark Function\n\nPress Any Key to Exit..." << endl;
		   getch();
		   exit(0);
	}

}

/* Fitness Function */
float Fitness(float a[])
{
	switch(BENCHMARK)
	{
		case 1: //No.1 - Sphere Function +-5.12
		   return Sphere(a);
		   
		case 2:	//No.2 - Ackley Function +-32.768
		   return Ackley(a);
		   
		case 3:	//No.3 - Rastrigin Function +-5.12
		   return Rastrigin(a);
		   
		case 4:	//No.4 - Zakharov Function -5, +10
		   return Zakharov(a);
		   
		case 5:	//No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
		   return AxisParallel(a);
		   
		case 6:	//No.6 - Griewank Function +-600
		   return Griewank(a);
		   
		case 7:	//No.7 - Sum of Different Powers function +-1.00
		   return SumOfDifferentPowers(a);
		   
		case 8: //No.8 - Rotated Hyper-Ellipsoid Function +-65.536
			return Rotated(a);
			
		case 9: //No.9 - Schwefel 2.22 Function +-5 [Updated]
		   return Schwefel(a);
		   
		case 10: //No.10 - Exponential function Function +-1.00 [f(x) = -1]
		   return Exponential(a);
		   
		default: //ERROR
		   cout << "Error: Invalid Benchmark Function\n\nPress Any Key to Exit..." << endl;
		   getch();
		   exit(0);
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
//Selection Function
//Done By: Yeap Chun Hong 2206352
//------------------------------------------------------------------------------------------------------------------------------
int RouletteWheelSelection2(float fitness[])
{
	float inverseFit[pSize];
	float totalFit = 0;
	
	for (int i = 0; i < pSize; i++)
	{
		inverseFit[i] = 1/fitness[i]; //inverse so that smaller fitness value will have a higher portion
		totalFit += inverseFit[i];
		//cout <<fitness[i] <<"\t" <<inverseFit[i] <<endl;
	}
	//cout <<totalFit <<endl;
	
	//calculate the cumulative probability and normalise it to [0,1]
	float cumulativeProbability[pSize];
    cumulativeProbability[0] = inverseFit[0] / totalFit;
	//cout <<"cumulativeProbability" << cumulativeProbability[0] <<endl;
    for(int i = 1; i < pSize; i++)
    {
        cumulativeProbability[i] = cumulativeProbability[i-1] + (inverseFit[i] / totalFit);
        //cout <<i << "\t"<<cumulativeProbability[i] <<endl;
    }
    
	// generate number from 0 to 1
	float spin = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);	//cout <<"spin: " <<spin <<endl;
	//cout << spin <<endl;
	
	for(int i = 0; i < pSize; i++)
    {
        if(spin <= cumulativeProbability[i])
        {
        	//cout << "chosen: "<<i<<endl;
            return i;
        }
    }

	// In case of rounding errors, return the last individual
	return pSize;
}

void RouletteWheelSelection(float fitness[], int& parent1, int& parent2)
{
	parent1 = RouletteWheelSelection2(fitness);
	parent2 = RouletteWheelSelection2(fitness);
	
	//cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;
    //getch();
}

int tournamentSize = 5;
int TournamentSelection2(float fitness[], int tournamentSize)
{
	int best = -1;
	float bestFitness = pow(999,30);
	bool selectedIndices[pSize] = {false};

//Randomly select unique individuals and perform tournament  
	for(int i = 0; i < tournamentSize; i++)
	{
		int index;		
		do 
		{
       		index = rand() % pSize;
    	} while (selectedIndices[index]);
		//cout <<"index "<< index <<endl;
	
    	selectedIndices[index] = true;

		//cout <<"fitness[index] "<< fitness[index]<<" bestFitness "<< bestFitness <<endl;
		if(fitness[index] < bestFitness)
		{
			best = index;
			bestFitness = fitness[index];
			//cout <<"best "<< best<<" bestFitness "<< bestFitness <<endl;
		}
	}
	//getch();
	return best;
}

void TournamentSelection(float fitness[], int tournamentSize, int& parent1, int& parent2)
{
	// Select first parent
    parent1 = TournamentSelection2(fitness, tournamentSize);

    // Select second parent
    parent2 = TournamentSelection2(fitness, tournamentSize);
    
    //getch();
}

int LinearRankingSelection2(float fitness[])
{
	//selection pressure
	float max = 1.1;
	float min = 2 - max; 
	
	
	//sort fitness value in descending order
	sort(fitness, fitness+ pSize, greater<float>());
	float probability[pSize];
	float totalFit = 0;
	for(int i = 0; i < pSize; i++)
    {
    	probability[i] = (min + (max - min) * i / (pSize -1) ) / pSize;
    	totalFit += probability[i];
        
    }
    
    //calculate the cumulative probability and normalise it to [0,1]
	float cumulativeProbability[pSize];
    cumulativeProbability[0] = probability[0] / totalFit;
    //cout <<"i" <<"\t"<< "Fitness"<<"\t"<< "Probability"<< "\t" << "Cumulative Prob" << endl;
    //cout <<"0" <<"\t"<< fitness[0]<<"\t"<< probability[0]<< "\t" <<cumulativeProbability[0] << endl;

	//cout <<"cumulativeProbability" << cumulativeProbability[0] <<endl;
    for(int i = 1; i < pSize; i++)
    {
        cumulativeProbability[i] = cumulativeProbability[i-1] + (probability[i] / totalFit);
        //cout <<i <<"\t"<< fitness[i]<<"\t"<< probability[i]<< "\t" <<cumulativeProbability[i] << endl;
        //cout <<i << "\t"<<cumulativeProbability[i] <<endl;
    }
    
	// generate number from 0 to 1
	float spin = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);	//cout <<"spin: " <<spin <<endl;
	//cout << spin <<endl;
	
	for(int i = 0; i < pSize; i++)
    {
        if(spin <= cumulativeProbability[i])
        {
        	//cout << "chosen: "<<i<<endl;
        	//getch();
            return i;
        }
    }

	// In case of rounding errors, return the last individual
	return pSize;
}

void LinearRankingSelection(float fitness[], int& parent1, int& parent2)
{
	// Select first parent
    parent1 = LinearRankingSelection2(fitness);

    // Select second parent
    parent2 = LinearRankingSelection2(fitness);
    
    //getch();
}


//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
//Crossover Function
//Done By: Loh Chia Heung 2301684
//------------------------------------------------------------------------------------------------------------------------------
void UniformCrossover(float chromosome[][dimension], float paroff[][dimension], float dcp, int p1, int p2) {
    float gcp = (rand() % 1000);
    gcp = gcp / 1000;
    
    if (gcp <= dcp) {
        // Perform Uniform Crossover
        for (int j = 0; j < dimension; j++) {
            if (rand() % 2 == 0) {
                // 0 --> Gene is copied from Parent 1
                paroff[2][j] = chromosome[p1][j];
                paroff[3][j] = chromosome[p2][j];
            } else {
                // 1 --> Gene is copied from Parent 2
                paroff[2][j] = chromosome[p2][j];
                paroff[3][j] = chromosome[p1][j];
            }
        }
    } else {
        // No crossover, directly copy parents to offspring
        for (int j = 0; j < dimension; j++) {
            paroff[2][j] = chromosome[p1][j]; // Offspring1 --> Parent 1 
            paroff[3][j] = chromosome[p2][j]; // Offspring2 --> Parent 2
        }
    }
}

void ShuffleCrossover(float chromosome[][dimension], float paroff[][dimension], float dcp, int p1, int p2) {
   
    float gcp = (rand() % 1000);
    gcp = gcp / 1000;
    
    if (gcp <= dcp) {
        float temp[dimension];

        // Shuffle chromosomes between parents
        for (int j = 0; j < dimension; j++) {
            temp[j] = (rand() % 2 == 0) ? chromosome[p1][j] : chromosome[p2][j];
        }

        // Assign shuffled chromosomes to offspring
        for (int j = 0; j < dimension; j++) {
            paroff[2][j] = temp[j]; 
            paroff[3][j] = (rand() % 2 == 0) ? chromosome[p1][j] : chromosome[p2][j]; 
        }
    } else {
        // No crossover --> offspring are exact copies of parents
        for (int j = 0; j < dimension; j++) {
            paroff[2][j] = chromosome[p1][j]; // Offspring1 --> Parent 1
            paroff[3][j] = chromosome[p2][j]; // Offspring2 --> Parent 2
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
//Mutation Function
//Done By: Brandon Ting En Junn 2101751
//------------------------------------------------------------------------------------------------------------------------------
void ReversingMutation()
{
	gmp = (rand() % 1000000);
   gmp = gmp / 1000000;
	
	//Parent 1
	if (gmp <= dmp)
	{
		mb1 = getrandom(0, dimension - 1);
		
		mb1v = paroff[2][mb1];
		
		if (mb1 == 0) {
			//Parent 1
			paroff[2][mb1] = paroff[2][dimension - 1];
			paroff[2][dimension - 1] = mb1v;
			
		} else {
			//Parent 1
			paroff[2][mb1] = paroff[2][mb1 - 1];
			paroff[2][mb1 - 1] = mb1v;
			
		}
	}
	
	gmp = (rand() % 1000000);
   gmp = gmp / 1000000;
	
	//Parent 2
	if (gmp <= dmp)
	{
		mb2 = getrandom(0, dimension - 1);
		
		mb2v = paroff[2][mb2];
		
		if (mb2 == 0) {
			//Parent 1
			paroff[3][mb2] = paroff[3][dimension - 1];
			paroff[3][dimension - 1] = mb2v;
			
		} else {
			//Parent 1
			paroff[3][mb1] = paroff[3][mb2 - 1];
			paroff[3][mb1 - 1] = mb2v;
			
		}
	}
}

void RandomMutation()
{
	//Parent 1 & Parent 2
	for (int i = 2; i < 4; i++)
	{
		gmp = (rand() % 1000000);
   	gmp = gmp / 1000000;
   
   	mb1 = getrandom(0 , dimension - 1);
   	mb2 = getrandom(0 , dimension - 1);
   	
   	if (gmp <= dmp) {
   		r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
   		paroff[i][mb1] = r;
		}
		
		gmp = (rand() % 1000000);
   	gmp = gmp / 1000000;
		
		if (gmp <= dmp) {
			r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
			paroff[i][mb2] = r;
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
//Replacement Function
//Done By: Ling Ji Xiang 2104584
//------------------------------------------------------------------------------------------------------------------------------
void WeakParentReplacement(float chromosome[pSize][dimension], float fit[pSize], float paroff[4][dimension], float tfit[4], int parent1, int parent2)
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

void BothParentReplacement(float chromosome[pSize][dimension], float fit[pSize], float paroff[4][dimension], float tfit[4], int parent1, int parent2)
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

//------------------------------------------------------------------------------------------------------------------------------
//External Function
//------------------------------------------------------------------------------------------------------------------------------
void resetExperiment()
{
	for (int i = 0; i < pSize; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
		 	chromosome[i][j] = 0;
		 		
		 	if (i < 4) {
		 		paroff[i][j] = 0;
		 		tfit[i] = 0;
			}
			
			if (i == 0) {
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
}

std::string to_string(int number) {
    if (number == 0) return "0";
    std::string result;
    bool isNegative = (number < 0);
    if (isNegative) number = -number;

    while (number > 0) {
        result += '0' + (number % 10);
        number /= 10;
    }

    if (isNegative) result += '-';
    std::reverse(result.begin(), result.end());
    return result;
}
//------------------------------------------------------------------------------------------------------------------------------

int main()
{
	srand(time(0));
	initialiseBenchmark_Range();
	
	// Directory e.g. GA1
	string gaDirectory = "GA" + to_string(BENCHMARK);
	CreateDirectory(gaDirectory.c_str(), NULL);
	
	// Directory e.g. GA1/GA1-Sphere
	string benchmarkDirectory = gaDirectory + "\\\\" + gaDirectory + benchmarkFunction;
	CreateDirectory(benchmarkDirectory.c_str(), NULL);
	
	//---------------------------------------------------------------------------------------------------------------------------
	//Experiment
	//---------------------------------------------------------------------------------------------------------------------------
	for (int i = 0; i < EXPERIMENT; i++)
	{
		// File e.g. GA1/GA1-Sphere/GA1-Sphere-Result1.txt
		string outfile1 = benchmarkDirectory + "\\\\" + gaDirectory + benchmarkFunction + "-Result" + to_string(i + 1) + ".txt";
		ofstream outfileo1(outfile1.c_str(),ios::trunc);
		
		cout << "Experiment " << i + 1 << "..."<< endl;
	
		//CPU Time
		clock_t start, end;
		start = clock();
//		srand(time(0));
    
		//---------------------------------------------------------------------------------------------------------------------------
		//Generate Population
		//---------------------------------------------------------------------------------------------------------------------------
		#if DEMO == 1
   		cout << "Generate Population" << endl;
		#endif
		
   		for(int i = 0 ; i < pSize ; i++)
   		{
      		for(int j = 0 ; j < dimension ; j++)
      		{
         		r = getrandom(-rangeMin,rangeMax);
         		r = r / rangeDiv;
         		chromosome[i][j] = r; 
      		}

			#if DEMO == 1
	      	if (i == 0 || i == pSize - 1) {
	      		cout<<"Chromosome "<<i+1<<endl;
		    	for(int j = 0 ; j < dimension ; j++)
		    	{
		        	cout<<setprecision(6)<<chromosome[i][j]<<"\t";
		      	}      
		      	cout<<endl<<endl;
			}
			#endif
   		}
//		getch();
		//---------------------------------------------------------------------------------------------------------------------------

		//---------------------------------------------------------------------------------------------------------------------------
		//Fitness Evaluation 1
		//---------------------------------------------------------------------------------------------------------------------------
		for(int i = 0 ; i < pSize ; i++)
		{
			fit[i] = Fitness(chromosome[i]);
//	      	cout<<setprecision(6)<<fit[i]<<endl;
			sumFit = 0;
		}
		
		lFv = pow(999,30);
   
		for(int i = 0; i < pSize; i++)
	    {
	       if(fit[i] < lFv)
	       {
	          lFv = fit[i];
	          lFvIndex = i;
	       }
	    }  
    
	    lowestGeneFV = lFv;
//	    cout<<"The best is chromosome "<<lFvIndex+1<<" with fitness of "<<lowestGeneFV<<endl;
	    for(int j = 0; j < dimension; j++)
	    {
	      lowestGene[j] = chromosome[lFvIndex][j]; 
	      //cout<<lowestGene[j]<<" ";
	    }  
//	    cout<<endl;
//	    getch();
//		outfileo1<<"Gen\tMinimum"<<endl;
		//---------------------------------------------------------------------------------------------------------------------------
		
		//---------------------------------------------------------------------------------------------------------------------------
		//Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
		//---------------------------------------------------------------------------------------------------------------------------
		for(int i = 0 ; i < gen ; i++)
		{
	   		//------------------------------------------------------------------------------------------------------------------------
      		//Selection Operation
      		//Done By: Yeap Chun Hong 2206352
      		//------------------------------------------------------------------------------------------------------------------------
			#if DEMO == 1
      		cout << "Selection Operation" << endl;
			#endif
			
			int parent1 = 0, parent2 = 0;
			
			
			//************************************************************************************************************************
     		/* Roulette Wheel Selection */
//    		RouletteWheelSelection(fit, parent1, parent2);
     		
     		/* Tournament Selection */
//			TournamentSelection(fit, tournamentSize, parent1, parent2);

			LinearRankingSelection(fit, parent1, parent2);
			
			//************************************************************************************************************************
//    cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;


			tfit[0] = fit[parent1]; 
			tfit[1] = fit[parent2]; 
      
			for(int j = 0 ; j < dimension ; j++)
			{
				paroff[0][j] = chromosome[parent1][j];
				paroff[1][j] = chromosome[parent2][j];
			}

			#if DEMO == 1
      		//Selected Parent 1
			cout<<"Chromosome "<<parent1 + 1<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
		   		cout<<setprecision(6)<<chromosome[parent1][j]<<"\t";
			}      
			cout<<endl<<endl;
		
			//Selected Parent 2
			cout<<"Chromosome "<<parent2 + 1<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
		   	cout<<setprecision(6)<<chromosome[parent2][j]<<"\t";
			}      
			cout<<endl<<endl;
			#endif
      		//------------------------------------------------------------------------------------------------------------------------
      
			//------------------------------------------------------------------------------------------------------------------------
			//Crossover Operation
			//Done By: Loh Chia Heung 2301684
			//------------------------------------------------------------------------------------------------------------------------
			#if DEMO == 1
      		cout << "Crossover Operation" << endl;
			#endif
			
			
			//************************************************************************************************************************
			/* Uniform Crossover */
			UniformCrossover(chromosome, paroff, dcp, parent1, parent2);
		
			/* Shuffle Crossover */
//			ShuffleCrossover(chromosome, paroff, dcp, parent1, parent2);
			//************************************************************************************************************************


			#if DEMO == 1
      		//Crossover Child 1
			cout<<"Child 1"<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
		   		cout<<setprecision(6)<<paroff[2][j]<<"\t";
			}      
			cout<<endl<<endl;
		
			//Crossover Child 2
			cout<<"Child 2"<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
			   cout<<setprecision(6)<<paroff[3][j]<<"\t";
			}      
			cout<<endl<<endl;
			#endif
      		//------------------------------------------------------------------------------------------------------------------------
      
      		//------------------------------------------------------------------------------------------------------------------------
      		//Mutation Operation
      		//Done By: Brandon Ting En Junn 2101751
      		//------------------------------------------------------------------------------------------------------------------------
			#if DEMO == 1
      		cout << "Mutation Operation" << endl;
			#endif
			
			
			//************************************************************************************************************************
			/* Reversing Mutation */
			ReversingMutation();
		
			/* Random Mutation */
//			RandomMutation();
			//************************************************************************************************************************


			#if DEMO == 1
      		//Mutation Child 1
			cout<<"Child 1"<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
			   cout<<setprecision(6)<<paroff[2][j]<<"\t";
			}      
			cout<<endl<<endl;
			
			//Mutation Child 2
			cout<<"Child 2"<<endl;
			for(int j = 0 ; j < dimension ; j++)
			{
			   cout<<setprecision(6)<<paroff[3][j]<<"\t";
			}      
			cout<<endl<<endl;
			#endif
      		//------------------------------------------------------------------------------------------------------------------------
      
			//------------------------------------------------------------------------------------------------------------------------
			//Fitness Evaluation 2
			//------------------------------------------------------------------------------------------------------------------------
			fit1 = Fitness(paroff[2]);
			sumFit = 0;
	       
			fit2 = Fitness(paroff[3]);
			sumFit = 0;
	
			tfit[2] = fit1; 
			tfit[3] = fit2;  
			//------------------------------------------------------------------------------------------------------------------------
      
      		//------------------------------------------------------------------------------------------------------------------------
      		//Replacement Operation
      		//Done By: Ling Ji Xiang 2104584
      		//------------------------------------------------------------------------------------------------------------------------
      		
      		
			//************************************************************************************************************************
			/* Weak Parent Replacement */
			WeakParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
		
			/* Both Parent Replacement */
//			BothParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
			//************************************************************************************************************************


      		//------------------------------------------------------------------------------------------------------------------------
      
      		lFv = pow(999,30);
			for(int j = 0; j < pSize; j++)
			{
				if(fit[j] < lFv)
	         	{
		            lFv = fit[j];
		            lFvIndex = j;
	         	}
	         
				//cout<<lFv<<" "<<lowestGeneFV<<endl;
				if(lFv < lowestGeneFV)
		        {
					lowestGeneFV = lFv;
					
					for(int k = 0; k < dimension; k++)
		            {
		            	lowestGene[k] = chromosome[lFvIndex][k]; 
		                //cout<<lowestGene[k]<<" ";
		            }  
		        }
	      	}  
	
	      	fit1 = 0;
	      	fit2 = 0;
	       
	      	outfileo1<<setprecision(6)<<lFv<<endl;
    	} //Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)] LOOP
    	//------------------------------------------------------------------------------------------------------------------------
    
	    lFv = pow(999,30);
	    for(int j = 0; j < pSize; j++)
	    {
	       if(fit[j] < lFv)
	       {
	          lFv = fit[j];
	          lFvIndex = j;
	       }
	    }  

		#if DEMO == 1
    	cout<<endl<<endl;
		#endif
    	outfileo1<<endl<<endl;
    	cout<<"Result"<<endl;
    	
    	//Output Final Generation Lowest
//    	cout<<setprecision(6)<<lFv<<" "<<lFvIndex+1<<endl<<endl;
//	    for(int j = 0 ; j < dimension ; j++)
//	    {
//	       cout<<setprecision(6)<<chromosome[lFvIndex][j]<<"\t";
//	       outfileo1<<setprecision(6)<<chromosome[lFvIndex][j]<<"\n";
//	    }
		
		//Output Lowest Gene in Experiment
    	cout<<setprecision(6)<<lowestGeneFV<<endl<<endl;
    	for(int j = 0 ; j < dimension ; j++)
    	{
    		cout<<setprecision(6)<<lowestGene[j]<<"\t";
       	outfileo1<<setprecision(6)<<lowestGene[j]<<"\n";
    	}
	    
	   cout<<endl;
	   outfileo1<<endl;
	   end = clock();
	   cout<<"Time required for execution: "<< (double)(end-start)/CLOCKS_PER_SEC<<" seconds."<<"\n\n";  
	   outfileo1<<(double)(end-start)/CLOCKS_PER_SEC<<"\n\n";
	   outfileo1.close();
		cout<<endl<<endl;
	 
		resetExperiment();
	} //Experiment LOOP
	//------------------------------------------------------------------------------------------------------------------------
	cout << "SUCCESSFULLY COMPLETED...\nPress Any Key to Exit..." << endl;
	getch();
	
	return 0;   
}
