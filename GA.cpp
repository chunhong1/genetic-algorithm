#include <iostream>
#include <conio.h>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <math.h>
#define getrandom(min,max) ((rand()%(((max)+1)-(min)))+(min))   
#define gen 2000           //number of iterations (number of generations)
#define pSize 40           //number of chromosomes (population size)
#define dimension 30       //number of bits (dimension size)
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
int rangeMin = 5120, rangeMax = 5120, rangeDiv = 1000;

int lFvIndex = 0;
float lFv = 0;

string outfile1 = ".\\GAResult1.txt"; 
ofstream outfileo1(outfile1.c_str(),ios::trunc);

//------------------------------------------------------------------------------------------------------------------------------
//Please write comment and tag your name while doing the coding
//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
//Fitness Function
//------------------------------------------------------------------------------------------------------------------------------
float Fitness(float a[])
{
   //---------------------------------------------------------------------------------------------------------------------------
   //Insert Benchmark Functions
   //---------------------------------------------------------------------------------------------------------------------------
   //Sphere +-5.12
   for(int j = 0 ; j < dimension ; j++) 
   {
      fv = pow(a[j],2);
      sumFit = sumFit + fv;
   }
   
   return sumFit;
}

int main()
{  
   //CPU Time
   clock_t start, end;
   start = clock();
   srand(time(0));
    
   //---------------------------------------------------------------------------------------------------------------------------
   //Generate Population
   //---------------------------------------------------------------------------------------------------------------------------
   for(int i = 0 ; i < pSize ; i++)
   {
      for(int j = 0 ; j < dimension ; j++)
      {
         r = getrandom(-rangeMin,rangeMax);
         r = r / rangeDiv;
         chromosome[i][j] = r; 
      }
      
      cout<<"Chromosome "<<i+1<<endl;
      for(int j = 0 ; j < dimension ; j++)
      {
         cout<<setprecision(6)<<chromosome[i][j]<<"\t";
      }      
      cout<<endl<<endl;
   }
   getch();
   
   //---------------------------------------------------------------------------------------------------------------------------
   //Fitness Evaluation
   //---------------------------------------------------------------------------------------------------------------------------
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
   for(int i = 0 ; i < gen ; i++)
   {
      //------------------------------------------------------------------------------------------------------------------------
      // Can modify starting from here (GA, DE, PSO)
      //------------------------------------------------------------------------------------------------------------------------

	  //------------------------------------------------------------------------------------------------------------------------
      //Selection Operation (Random Selection)
      //------------------------------------------------------------------------------------------------------------------------
      int parent1 = rand() % pSize;
      
      redo:
      int parent2 = rand() % pSize;
      
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
      //Crossover Operation (Single Point Crossover)
      //------------------------------------------------------------------------------------------------------------------------
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
      //Mutation Operation (Interchanging Mutation)
      //------------------------------------------------------------------------------------------------------------------------
      gmp = (rand() % 1000000);
      gmp = gmp / 1000000;
      mb1 = getrandom(0 , dimension-1);
      
      redo2:
      mb2 = getrandom(0 , dimension-1);
      
      if(mb2 == mb1)
      {
         goto redo2;           
      }    
      
      if(gmp<=dmp)
      {
         mb1v = paroff[2][mb1];
         mb2v = paroff[2][mb2];  
         paroff[2][mb1] = mb2v;
         paroff[2][mb2] = mb1v;
         mb1v = paroff[3][mb1];
         mb2v = paroff[3][mb2];  
         paroff[3][mb1] = mb2v;
         paroff[3][mb2] = mb1v;
      }
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
      //Replacement Operation (Random Replacement)  
      //------------------------------------------------------------------------------------------------------------------------
      rp1 = getrandom(0,3);
      
      redo3:
      rp2 = getrandom(0,3);
      
      if(rp2 == rp1)
      {
         goto redo3;           
      }    
    
      for(int j = 0 ; j < dimension ; j++)
      {
         chromosome[parent1][j] =  paroff[rp1][j];  
         chromosome[parent2][j] =  paroff[rp2][j];            
      }
      
      fit[parent1] = tfit[rp1];
      fit[parent2] = tfit[rp2];
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

    cout<<endl<<endl;
    outfileo1<<endl<<endl;
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
