#include <iostream>
#include <conio.h>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <math.h>
using namespace std;

#define _USE_MATH_DEFINES
const float pi = M_PI;

/* EXPERIMENT */
#define getrandom(min,max) ((rand()%(((max)+1)-(min)))+(min)) //random integer between min and max (inclusive)
#define pSize 40           											  //population size 		  (no. of chromosomes)
#define dimension 30       											  //dimenstion size		  (number of bits)
#define gen 2000           											  //number of generations (no. of iterations)

/* CHROMOSOME */
float chromosome[pSize][dimension]; 								  //chromosome
float paroff[4][dimension];         								  //parent and offspring (child) -> parent1[0], parent2[1], offspring1[2], offspring2[3]
float fit[pSize];                   								  //fitness value for each chromosome

float r = 0;																  //r = random
// Crossover
float dcp = 0.7, gcp = 0;												  //dcp = default crossover probability, gcp = generated crossover probability
int crb = 0;																  //crb = crossover breakpoint
// Mutation
float dmp = 0.01, gmp = 0, gmp1 = 0, gmp2 = 0;					  //dmp = default mutation probability, gmp = generated mutation probability
int mb1 = 0, mb2 = 0;													  //mb = mutation positions
float mb1v = 0 , mb2v = 0;												  //mbv = mutation values
// Replacement
int rp1 = 0, rp2 = 0;													  //rp = replacement indices

float fv = 0, sumFit=0;													  //fv = fitness value
float fit1 = 0 , fit2 = 0;
float tfit[4];

int lFvIndex = 0;
float lFv = 0; 															  //lFv = lowest fitness value

/* RANGE */
int rangeMin = 5120, rangeMax = 5120, rangeDiv = 1000;		  //r = getrandom(-rangeMin,rangeMax) / rangeDiv

/* TEXT FILE */
string outfile1 = ".\\GAResult1.txt"; 
ofstream outfileo1(outfile1.c_str(),ios::trunc);

int tournamentSize = 5;

//------------------------------------------------------------------------------------------------------------------------------
//Fitness Function
//------------------------------------------------------------------------------------------------------------------------------
//10 Benchmark Functions
float Sphere(float a[]);
float Ackley(float a[]);
float Rastrigin(float a[]);
float Zakharov(float a[]);
float AxisParallel(float a[]);

float Griewank(float a[]);
float SumOfDifferentPowers(float a[]);
float Rotated(float a[]);
float Schwefel(float a[]);
float Exponential(float a[]);

/* Fitness */
float Fitness(float a[])
{
   //No.1 - Sphere  Function +-5.12
   return Sphere(a);
}

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
//Written By: Yeap Chun Hong 2206352
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
//Written By: Yeap Chun Hong 2206352
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
//Written By: Brandon Ting En Junn
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
//Written By: Loh Chia Heung 2301684
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
//Written By: Loh Chia Heung 2301684
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
//Written By: Yeap Chun Hong 2206352
float SumOfDifferentPowers(float a[])
{
   for(int j = 0 ; j < dimension ; j++) 
   {
   	fv = pow(fabs(a[j]),j+1);
      sumFit = sumFit + fv;
   }
   return sumFit;
}

//No.8 - Rotated Hyper-Ellipsoid Function +-65.536
//Written By: Brandon Ting En Junn 2101751
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

//No.9 - Schewefel 2.22 Function +-100 [ERROR]
//Written By: Ling Ji Xiang 2104584
float Schwefel(float a[])
{
	for(int i = 0; i < dimension; i++)
	{
		float sum = 0.0, product = 1.0;
		float absolute = abs(a[i]);
		sum += absolute;
		product *= absolute;
	}
	return sum + product;
}

//No.10 - Exponential function Function +-1.00 [REVISE]
//Written By: Ling Ji Xiang 2104584
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

//------------------------------------------------------------------------------------------------------------------------------
//Operation Technique
//------------------------------------------------------------------------------------------------------------------------------
/*/////////////////////////////////////////
	Selection Function
	Written By: Yeap Chun Hong 2206352
*//////////////////////////////////////////
int RouletteWheelSelection(float fitness[])
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

int TournamentSelection( float fitness[], int tournamentSize)
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


/*/////////////////////////////////////////
	Crossover Function
	Written By: Loh Chia Heung 2301684
*//////////////////////////////////////////



/*/////////////////////////////////////////
	Mutation Function
	Written By: Brandon Ting En Junn 2101751
*//////////////////////////////////////////
void ReversingMutation()
{
	gmp1 = (rand() % 1000000);
   gmp1 = gmp1 / 1000000;
   
   gmp2 = (rand() % 1000000);
   gmp2 = gmp2 / 1000000;
	
	//Parent 1
	if (gmp1 <= dmp)
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
	
	//Parent 2
	if (gmp2 <= dmp)
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
		gmp1 = (rand() % 1000000);
   	gmp1 = gmp1 / 1000000;
   
   	gmp2 = (rand() % 1000000);
   	gmp2 = gmp2 / 1000000;
   
   	mb1 = getrandom(0 , dimension - 1);
   	mb2 = getrandom(0 , dimension - 1);
   	
   	if (gmp1 <= dmp) {
   		r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
   		paroff[i][mb1] = r;
		}
		
		if (gmp2 <= dmp) {
			r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
			paroff[i][mb2] = r;
		}
	}
}


/*/////////////////////////////////////////
	Replacement Function
	Written By: Ling Ji Xiang 2104584
*//////////////////////////////////////////
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


int main()
{  
	cout<<"GA001"<<endl;
   //CPU Time
   clock_t start, end;
   start = clock();
   srand(time(0));
    
   //---------------------------------------------------------------------------------------------------------------------------
   //Generate Population
   //---------------------------------------------------------------------------------------------------------------------------
   for(int i = 0 ; i < pSize ; i++)
   {
   	//Generate all random gene values (j) for chromosome (i)
      for(int j = 0 ; j < dimension ; j++)
      {
         r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
         chromosome[i][j] = r; 
      }
      
      //Display chromosome (i) with all gene values, d.p. = 6 (j)
      cout<<"Chromosome "<<i+1<<endl;
      for(int j = 0 ; j < dimension ; j++)
      {
         cout<<setprecision(6)<<chromosome[i][j]<<"\t";
      }
//      cout<<setprecision(6)<<Fitness(chromosome[i])<<endl;
//      sumFit = 0;
      cout<<endl<<endl;
   }
//   getch();
   
   //---------------------------------------------------------------------------------------------------------------------------
   //Fitness Evaluation
   //---------------------------------------------------------------------------------------------------------------------------
	cout<<"Evaluating Fitness..."<<endl;
	float total_sumFit = 0;
	for(int i = 0 ; i < pSize ; i++)
   {
      fit[i] = Fitness(chromosome[i]);
      cout<<setprecision(6)<<fit[i]<<endl;
      sumFit = 0;
   }
  
   //outfileo1<<"Gen\tMinimum"<<endl;
   
   
   //---------------------------------------------------------------------------------------------------------------------------
   //Termination Criteria (Maximum Generation)
   //---------------------------------------------------------------------------------------------------------------------------
   cout<<endl<<endl<<"Searching "<<gen<<" Generations...";
	for(int i = 0 ; i < gen ; i++)
   {
      //------------------------------------------------------------------------------------------------------------------------
      // Can modify starting from here (GA, DE, PSO)
      //------------------------------------------------------------------------------------------------------------------------


	   //------------------------------------------------------------------------------------------------------------------------
      //Selection Operation (Roulette Wheel Selection / Tournament Selection)
      //Done By: Yeap Chun Hong 2206352
      //------------------------------------------------------------------------------------------------------------------------
      //int parent1 = rand() % pSize;
      //int parent1 = RouletteWheelSelection(fit);
      int parent1 = TournamentSelection(fit, tournamentSize);
      redo:
      //int parent2 = RouletteWheelSelection(fit);
      int parent2 = TournamentSelection(fit, tournamentSize);      
      
      if(parent2 == parent1)
      {
         goto redo;           
      }           
  
      tfit[0] = fit[parent1]; 
      tfit[1] = fit[parent2]; 
      
      for(int j = 0 ; j < dimension ; j++)
      {
         paroff[0][j] = chromosome[parent1][j];
         paroff[1][j] = chromosome[parent2][j];
      }
		//------------------------------------------------------------------------------------------------------------------------ 
      
      
      //------------------------------------------------------------------------------------------------------------------------
      //Crossover Operation (Uniform Crossover / Three Parent Crossover)
      //Done By: Loh Chia Heung 2301684
      //------------------------------------------------------------------------------------------------------------------------
      
      //Sample - Single Point Crossover
      gcp = (rand() % 1000);
      gcp = gcp / 1000;
      crb = getrandom(1 , dimension-2);
      
      if(gcp<=dcp)
      {
         for(int j = 0 ; j < crb; j++)
         {       
            paroff[2][j] = paroff[0][j];
            paroff[3][j] = paroff[1][j]; 
         } 
         for(int j = crb ; j < dimension; j++)
         {       
            paroff[2][j] = paroff[1][j];
            paroff[3][j] = paroff[0][j]; 
         }      
      }  
      else
      {
         for(int j = 0 ; j < dimension; j++)
         {       
            paroff[2][j] = paroff[0][j];
            paroff[3][j] = paroff[1][j]; 
         }     
      }
      //------------------------------------------------------------------------------------------------------------------------
      
      
      //------------------------------------------------------------------------------------------------------------------------
      //Mutation Operation (Reversing Mutation / Random Mutation)
      //Done By: Brandon Ting En Junn
      //------------------------------------------------------------------------------------------------------------------------
      //Sample - Interchanging Mutation
//		gmp = (rand() % 1000000);
//      gmp = gmp / 1000000;
//      mb1 = getrandom(0 , dimension-1);
//      
//      redo2:
//      mb2 = getrandom(0 , dimension-1);
//      
//      if(mb2 == mb1)
//      {
//         goto redo2;           
//      }    
//      
//      if(gmp<=dmp)
//      {
//         mb1v = paroff[2][mb1];
//         mb2v = paroff[2][mb2];  
//         paroff[2][mb1] = mb2v;
//         paroff[2][mb2] = mb1v;
//         mb1v = paroff[3][mb1];
//         mb2v = paroff[3][mb2];  
//         paroff[3][mb1] = mb2v;
//         paroff[3][mb2] = mb1v;
//      }
//		ReversingMutation();
		RandomMutation();
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
      //Replacement Operation
      //Done By: Ling Ji Xiang 2104584
      //------------------------------------------------------------------------------------------------------------------------
      // Apply the replacement strategies
      // Sample - WorstParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
      BothParentReplacement(chromosome, fit, paroff, tfit, parent1, parent2);
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
       
      outfileo1<<setprecision(6)<<lFv<<endl;
    }  
    
    /* RESULT */
    lFv = pow(999,30);
    for(int j = 0; j < pSize; j++)
    {
       if(fit[j] < lFv)
       {
          lFv = fit[j];
          lFvIndex = j;
       }
    }  

    cout<<endl;
    outfileo1<<endl<<endl;
    cout<<"Best Chromosome"<<endl;
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
    getch();
    return 0;   
}
