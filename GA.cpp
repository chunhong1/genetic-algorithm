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

//------------------------------------------------------------------------------------------------------------------------------
//Parameter Settings
//1	/* Benchmark Function */ 			- 1 to 10
//1.1 /* Range */
const int BENCHMARK = 1;

//2.	/* Selection */ 						- Roulette Wheel / Tournament
//2.1 /* Selection Operation */ 			- Roulette Wheel

//3.	/* Crossover */ 						- Uniform / Shuffle
//3.1 /* Crossover Operation */ 			- Uniform

//4.	/* Mutation */ 						- Reversing / Random
//4.1 /* Mutation Operation */ 			- Reversing

//5.	/* Replacement */ 					- Weak Parent / Both Parent
//5.1 /* Replacement Operation */ 		- Weak Parent

// /* Demo */									- 1 Generation (Iteration) + CMD Output
#define DEMO 0

// /* Experiment */							- Run 10 Times
#define EXPERIMENT 10
//------------------------------------------------------------------------------------------------------------------------------

#define getrandom(min,max) (static_cast<long long>(rand()) * (max - min + 1) / RAND_MAX) + min
#if DEMO == 1
	#define gen 1
#else
	#define gen 2000           //number of iterations (number of generations)
#endif
#define pSize 40           	//number of chromosomes (population size)
#define dimension 30       	//number of bits (dimension size)
using namespace std;

float chromosome[pSize][dimension]; //chromosome
float paroff[4][dimension];         //parent and offspring
float fit[pSize];                   //fitness value for each chromosome
float r = 0, dcp = 0.7, dmp = 0.01, gcp = 0, gmp = 0;
int crb = 0 , mb1 = 0 , mb2 = 0;
int rp1 = 0, rp2 = 0;
float mb1v = 0 , mb2v = 0;
float fv = 0, sumFit=0;
float fit1 = 0 , fit2 = 0;
float tfit[4];

int lFvIndex = 0;
float lFv = 0;

const float pi = M_PI;

int tournamentSize = 5;

//------------------------------------------------------------------------------------------------------------------------------
//Fitness Function
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
//Done By: Brandon Ting En Junn
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
    float sumFit = 0.0;
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

/* Range */
int rangeMin = 0, rangeMax = 0, rangeDiv = 1000;
void initialiseRange()
{
	switch(BENCHMARK)
	{
		case 1: //No.1 - Sphere Function +-5.12
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 2:	//No.2 - Ackley Function +-32.768
		   rangeMin = rangeMax = 32768;
		   break;
		   
		case 3:	//No.3 - Rastrigin Function +-5.12
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 4:	//No.4 - Zakharov Function -5, +10
		   rangeMin = 5000; rangeMax = 10000;
		   break;
		   
		case 5:	//No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
		   rangeMin = rangeMax = 5120;
		   break;
		   
		case 6:	//No.6 - Griewank Function +-600
		   rangeMin = rangeMax = 600000;
		   break;
		   
		case 7:	//No.7 - Sum of Different Powers function +-1.00
		   rangeMin = rangeMax = 1000;
		   break;
		   
		case 8: 	//No.8 - Rotated Hyper-Ellipsoid Function +-65.536
			rangeMin = rangeMax = 65536;
		   break;
		   
		case 9: 	//No.9 - Schwefel 2.22 Function +-5 [Updated]
		   rangeMin = rangeMax = 5000;
		   break;
		   
		case 10: //No.10 - Exponential function Function +-1.00 [f(x) = -1]
		   rangeMin = rangeMax = 1000;
		   break;
		   
		default: //ERROR
		   cout << "BENCHMARK FUNCTION ERROR" << endl;
		   getch();
		   exit(0);
	}

}


/* Benchmark Function */
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
		   
		case 8: 	//No.8 - Rotated Hyper-Ellipsoid Function +-65.536
			return Rotated(a);
			
		case 9: 	//No.9 - Schwefel 2.22 Function +-5 [Updated]
		   return Schwefel(a);
		   
		case 10: //No.10 - Exponential function Function +-1.00
		   return Exponential(a);
		   
		default: //ERROR
		   cout << "BENCHMARK FUNCTION ERROR" << endl;
		   getch();
		   exit(0);
	}
}

//------------------------------------------------------------------------------------------------------------------------------
///* Selection */
//Done By: Yeap Chun Hong 2206352
//------------------------------------------------------------------------------------------------------------------------------
int RouletteWheelSelection2(float fitness[])
{
	//calculate the total fitness
	float total_sumFit = 0;
	for (int i=0; i<pSize;i++)
	{
		total_sumFit += fitness[i];
	}
	
	//normalise fitness value
  float normFit[pSize];
  for(int i = 0 ; i < pSize ; i++)
  {
    normFit[i] = fit[i]/total_sumFit;  
  }
   
   // Calculate cumulative probabilities
  float cumulativeProb[pSize];
  cumulativeProb[0] = normFit[0];
    
  for (int i = 1; i < pSize; i++) 
	{
		cumulativeProb[i] = cumulativeProb[i - 1] + normFit[i];
  }
    
  // generate number 0 to 1
  float spin = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        
  // Select the individual based on the spin
  for (int i = 0; i < pSize; i++) 
	{
		if (spin < cumulativeProb[i]) 
		{
			return i;
    }
  }
  // In case of rounding errors, return the last individual
  return pSize - 1;
}

void RouletteWheelSelection(float fitness[], int& parent1, int& parent2)
{
	parent1 = RouletteWheelSelection2(fitness);
	parent2 = RouletteWheelSelection2(fitness);
	
	//cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;
    //getch();
}

int TournamentSelection2(float fitness[], int tournamentSize)
{
	int best = -1;
	float bestFitness = pow(999,-30);
	bool selectedIndices[pSize] = {false};

//Randomly select unique individuals and perform tournament  
	for(int i = 0; i < tournamentSize; i++)
	{
		int index;		
		do 
		{
       index = rand() % pSize;
    } while (selectedIndices[index]);

    selectedIndices[index] = true;

		if(fitness[index] > bestFitness)
		{
			best = index;
			bestFitness = fitness[index];
		}
	}
	
	return best;
}

void TournamentSelection(float fitness[], int tournamentSize, int& parent1, int& parent2)
{
	// Select first parent
    parent1 = TournamentSelection2(fitness, tournamentSize);

    // Select second parent
    parent2 = TournamentSelection2(fitness, tournamentSize);
    
    //cout <<"Parent1: "<<parent1<<"parent2 "<<parent2<<endl;
    //getch();
}

//------------------------------------------------------------------------------------------------------------------------------
///* Crossover */
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
///* Mutation */
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
///* Replacement */
//Done By: Ling Ji Xiang 2104584
//------------------------------------------------------------------------------------------------------------------------------
void WeakParentReplacement(float chromosome[pSize][dimension], float fit[pSize], float paroff[4][dimension], float tfit[4], int parent1, int parent2)
{
    if (fit[parent1] > fit[parent2]) //Check which parent has worse fitness value (higher fitness value)*
    {
        if (tfit[2] < fit[parent1]) //Check if fitness value of offspring1 is lower than parent1
        {
            for (int i = 0; i < dimension; i++)
            {
		//Replace genes of parent1 with offspring1
                chromosome[parent1][i] = paroff[2][i]; 
            }
            fit[parent1] = tfit[2]; //Update fitness of parent1
        }
        if (tfit[3] < fit[parent2]) //Check if fitness value of offspring2 is lower than parent2
        {
            for (int i = 0; i < dimension; i++)
            {
		//Replace genes of parent2 with offspring2
                chromosome[parent2][i] = paroff[3][i]; 
            }
            fit[parent2] = tfit[3]; //Update fitness of parent2
        }
    }
    else
    {
        if (tfit[2] < fit[parent2]) //Check if fitness value of offspring1 is lower than parent2
        {
            for (int i = 0; i < dimension; i++)
            {
		//Replace genes of parent2 with offspring1
                chromosome[parent2][i] = paroff[2][i]; 
            }
            fit[parent2] = tfit[2]; // Update fitness of parent2
        }
        if (tfit[3] < fit[parent1]) //Check if fitness value of offspring2 is lower than parent1
        {
            for (int i = 0; i < dimension; i++)
            {
		//Replace genes of parent1 with offspring2
                chromosome[parent1][i] = paroff[3][i]; 
            }
            fit[parent1] = tfit[3]; // Update fitness of parent1
        }
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

// Int to String Function
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

int main()
{
	srand(time(0));
	initialiseRange();
	
	//BENCHMARK LOOP
	
	//EXPERIMENT LOOP
	for (int i = 0; i < EXPERIMENT; i++)
	{
		cout << "Experiment " << i + 1 << "..."<< endl;
		string outfile1 = ".\\GAResult" + to_string(i + 1) + ".txt";
		ofstream outfileo1(outfile1.c_str(),ios::trunc);

   //CPU Time
   clock_t start, end;
   start = clock();
//	srand(time(0));
    
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
//   getch();
   
   //---------------------------------------------------------------------------------------------------------------------------
   //Fitness Evaluation
   //---------------------------------------------------------------------------------------------------------------------------
   for(int i = 0 ; i < pSize ; i++)
   {
      fit[i] = Fitness(chromosome[i]);
//      cout<<setprecision(6)<<fit[i]<<endl;
      sumFit = 0;
   }
  
//   outfileo1<<"Gen\tMinimum"<<endl;
   //---------------------------------------------------------------------------------------------------------------------------
   //Termination Criteria (Maximum Generation)
   //---------------------------------------------------------------------------------------------------------------------------
   for(int i = 0 ; i < gen ; i++)
   {
      //------------------------------------------------------------------------------------------------------------------------
      // Can modify starting from here (GA, DE, PSO)
      //------------------------------------------------------------------------------------------------------------------------

	   //------------------------------------------------------------------------------------------------------------------------
      ///* Selection Operation */
      //Done By: Yeap Chun Hong 2206352
      //------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
      cout << "Selection Operation" << endl;
#endif
		int parent1 = 0, parent2 = 0;
		
     	//Roulette Wheel Selection
     	RouletteWheelSelection(fit, parent1, parent2);
     
     	//Tournament Selection
//     	TournamentSelection(fit, tournamentSize, parent1, parent2);

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
      ///* Crossover Operation */
      //Done By: Loh Chia Heung 2301684
      //------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
      cout << "Crossover Operation" << endl;
#endif
		//Uniform Crossover
		UniformCrossover(chromosome, paroff, dcp, parent1, parent2);
		
		//Shuffle Crossover
//		ShuffleCrossover(chromosome, paroff, dcp, parent1, parent2);

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
      ///* Mutation Operation */
      //Done By: Brandon Ting En Junn 2101751
      //------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
      cout << "Mutation Operation" << endl;
#endif
		//Reversing Mutation
		ReversingMutation();
		
		//Random Mutation
//		RandomMutation();

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
      //Fitness Evaluation
      //------------------------------------------------------------------------------------------------------------------------
      fit1 = Fitness(paroff[2]);
      sumFit = 0;
       
      fit2 = Fitness(paroff[3]);
      sumFit = 0;

      tfit[2] = fit1; 
      tfit[3] = fit2;  
      //------------------------------------------------------------------------------------------------------------------------
      
      //------------------------------------------------------------------------------------------------------------------------
      ///* Replacement Operation */
      //Done By: Ling Ji Xiang 2104584
      //------------------------------------------------------------------------------------------------------------------------
		//Weak Parent Replacement
		WeakParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
		
		//Both Parent Replacement
//		BothParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
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

      fit1 = 0;
      fit2 = 0;
       
      outfileo1<<setprecision(6)<<lFv<<endl;;
    }  
    
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
    cout<<setprecision(6)<<lFv<<" "<<lFvIndex+1<<endl<<endl;
    
    for(int j = 0 ; j < dimension ; j++)
    {
       cout<<setprecision(6)<<chromosome[lFvIndex][j]<<"\t";
       outfileo1<<setprecision(6)<<chromosome[lFvIndex][j]<<"\n";
    }
    
    cout<<endl;
    outfileo1<<endl;
    end = clock();
    cout<<"Time required for execution: "<< (double)(end-start)/CLOCKS_PER_SEC<<" seconds."<<"\n\n";  
    outfileo1<<(double)(end-start)/CLOCKS_PER_SEC<<"\n\n";
    outfileo1.close();
	 cout<<endl<<endl;
	 
	} //EXPERIMENT LOOP
	 cout << "SUCCESSFULLY COMPLETED..." << endl;
    getch();
    return 0;   
}
