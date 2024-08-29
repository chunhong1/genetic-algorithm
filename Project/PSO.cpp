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
#include <windows.h>
using namespace std;

/******************************************************** CODE GUIDELINE ********************************************************
* Parameter Settings (PSO with constant w and vc)

- Benchmark Function
- External Function
* Update Velocity Function
* Update Position Function

-> Benchmark Automation
	-> Experiment Automation
		-> Generate Population
		-> Fitness Evaluation 1
		-> Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
			*> Update Velocity
			*> Update Position
			-> Fitness Evaluation 2
******************************************************** CODE GUIDELINE ********************************************************/

//------------------------------------------------------------------------------------------------------------------------------
// Parameter Settings (PSO with constant w and vc)
//------------------------------------------------------------------------------------------------------------------------------
#define DEMO 0

const double w = 0.7; 	   							//Inertia Weight
const double c1 = 2, c2 = 2;						//C Constant
//------------------------------------------------------------------------------------------------------------------------------
#define EXPERIMENT 10								//No. of Experiments

// #define getrandom(min,max) ((rand()%(((max)+1)-(min)))+(min))
#define getrandom(min, max) (((double)rand() / RAND_MAX) * ((max) - (min)) + (min))
#define gen 2000           							//number of iterations (number of generations)
#define pSize 40           							//number of particle (population size)
#define dimension 30       							//number of bits (dimension size)

double particlePosition[pSize][dimension] = {0}; 	//particle position
double particleVelocity[pSize][dimension] = {0};	//particle velocity
double particlePBest[pSize][dimension] = {0};       //particle PBest
double particlePBestFV[pSize] = {0};
double particleGBest[dimension] = {0};              //particle GBest
double particleGBestFV = 0;

double fit[pSize] = {0};                            //fitness value for each particle generation
double fv = 0, sumFit = 0;
double r = 0;
int GBestIndex = 0;									//GBest index of each particle generation for potential GBest

int BENCHMARK = 1;
string benchmarkFunction = "";
int rangeMin = 0, rangeMax = 0, rangeDiv = 1000;

double vcMin = 0, vcMax = 0; 	   					// Velocity Clamping (Velocity Limits)
double posMin = 0, posMax = 0;						// Position Limits

//------------------------------------------------------------------------------------------------------------------------------
// Benchmark Function
//------------------------------------------------------------------------------------------------------------------------------
// No.1 - Sphere Function +-5.12
double Sphere(double a[])
{
	sumFit = 0;

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
	sumFit = 0;
	double sum1 = 0, sum2 = 0;

	for (int j = 0; j < dimension; j++)
	{
		sum1 += pow(a[j], 2);
		sum2 += cos(2 * M_PI * a[j]);
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
	sumFit = 0;

	for (int j = 0; j < dimension; j++)
	{
		fv = (pow(a[j], 2)) - (10 * cos(2 * M_PI * a[j]));
		sumFit = sumFit + fv;
	}

	sumFit += 10 * dimension;

	return sumFit;
}

// No.4 - Zakharov Function -5, +10
// Done By: Brandon Ting En Junn 2101751
double Zakharov(double a[])
{
	sumFit = 0;
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
	sumFit = 0;

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
	sumFit = 0;
	double product = 1.0;

	for (int j = 0; j < dimension; j++)
	{
		sumFit += (a[j] * a[j]) / 4000.0;
		product *= cos(a[j] / sqrt(j + 1));
	}

	sumFit = sumFit - product + 1;

	return sumFit;
}

// No.7 - Sum of Different Powers function +-1.00
// Done By: Yeap Chun Hong 2206352
double SumOfDifferentPowers(double a[])
{
	sumFit = 0;

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
	sumFit = 0;

	for (int i = 0; i < dimension; i++)
	{
		fv = 0;

		for (int j = 0; j <= i; j++)
		{
			fv += pow(a[j], 2);
		}

		sumFit = sumFit + fv;
	}

	return sumFit;
}

// No.9 - Schwefel 2.22 Function +-5 [Updated]
// Done By: Ling Ji Xiang 2104584
double Schwefel(double a[])
{
	sumFit = 0;
	double absolute = 0, sum = 0, product = 1.0;

	for (int i = 0; i < dimension; i++)
	{
		absolute = fabs(a[i]);
		sum += absolute;

		product *= absolute;
	}

	sumFit = sum + product;

	return sumFit;
}

// No.10 - Exponential function Function +-1.00 [f(x) = -1]
// Done By: Ling Ji Xiang 2104584
double Exponential(double a[])
{
	sumFit = 0;

	// calculate the sum of squares
	for (int j = 0; j < dimension; j++)
	{
		fv = pow(a[j], 2);
		sumFit = sumFit + fv;
	}

	// Calculate the exponential of sum of squares
	sumFit = -exp(-0.5 * sumFit);

	return sumFit;
}

/* Benchmark_Range */
void initialiseBenchmark_Range()
{
	switch (BENCHMARK)
	{
	case 1: // No.1 - Sphere Function +-5.12
		benchmarkFunction = "-Sphere";
		rangeMin = rangeMax = 5120;

		vcMax = 5120 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 5120 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 2: // No.2 - Ackley Function +-32.768
		benchmarkFunction = "-Ackley";
		rangeMin = rangeMax = 32768;

		vcMax = 32768 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 32768 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 3: // No.3 - Rastrigin Function +-5.12
		benchmarkFunction = "-Rastrigin";
		rangeMin = rangeMax = 5120;

		vcMax = 5120 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 5120 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 4: // No.4 - Zakharov Function -5, +10
		benchmarkFunction = "-Zakharov";
		rangeMin = 5000;
		rangeMax = 10000;

		vcMax = 10000 / static_cast<double>(rangeDiv) / 2;
		vcMin = -(5000 / static_cast<double>(rangeDiv) / 2);

		posMax = 10000 / static_cast<double>(rangeDiv);
		posMin = -(5000 / static_cast<double>(rangeDiv));
		break;

	case 5: // No.5 - Axis Parallel Hyper-Ellipsoid Function +-5.12
		benchmarkFunction = "-AxisParallel";
		rangeMin = rangeMax = 5120;

		vcMax = 5120 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 5120 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 6: // No.6 - Griewank Function +-600
		benchmarkFunction = "-Griewank";
		rangeMin = rangeMax = 600000;

		vcMax = 600000 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 600000 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 7: // No.7 - Sum of Different Powers function +-1.00
		benchmarkFunction = "-SumOfDifferentPowers";
		rangeMin = rangeMax = 1000;

		vcMax = 1000 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 1000 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 8: // No.8 - Rotated Hyper-Ellipsoid Function +-65.536
		benchmarkFunction = "-Rotated";
		rangeMin = rangeMax = 65536;

		vcMax = 65536 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 65536 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 9: // No.9 - Schwefel 2.22 Function +-5 [Updated]
		benchmarkFunction = "-Schwefel";
		rangeMin = rangeMax = 5000;

		vcMax = 5000 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 5000 / static_cast<double>(rangeDiv);
		posMin = -posMax;
		break;

	case 10: // No.10 - Exponential function Function +-1.00 [f(x) = -1]
		benchmarkFunction = "-Exponential";
		rangeMin = rangeMax = 1000;

		vcMax = 1000 / static_cast<double>(rangeDiv) / 2;
		vcMin = -vcMax;

		posMax = 1000 / static_cast<double>(rangeDiv);
		posMin = -posMax;
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

void resetExperiment()
{
	for (int i = 0; i < pSize; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			particlePosition[i][j] = 0;
			particleVelocity[i][j] = 0;
			particlePBest[i][j] = 0;

			if (i == 0)
			{
				particleGBest[j] = 0;
			}
		}

		particlePBestFV[pSize] = 0;

		fit[i] = 0;
	}
	particleGBestFV = 0;

	fv = sumFit = 0;
	r = 0;
	GBestIndex = 0;
}
//------------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------------
// Update Velocity Function -> Refer w | vcMin, vcMax as lower and upper bounds | getrandom(0, 1000) / 1000;
// Done By: Brandon Ting En Junn 2101751
//------------------------------------------------------------------------------------------------------------------------------
double PSO_Constant_W_VC(double position, double velocity, double PBest, double GBest)
{
	// Constant W
	double v = w * velocity + c1 * static_cast<double>(getrandom(0, 1000)) / 1000 * (PBest - position) + c2 * static_cast<double>(getrandom(0, 1000)) / 1000 * (GBest - position);

	// Velocity Clamping
	if (v < vcMin) {
		return vcMin;
	} else if (v > vcMax) {
		return vcMax;
	} else {
		return v;
	}
}

void UpdateVelocity()
{
	for (int i = 0; i < pSize; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			particleVelocity[i][j] = PSO_Constant_W_VC(particlePosition[i][j], particleVelocity[i][j], particlePBest[i][j], particleGBest[j]);
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------
// Update Position Function -> Refer posMin, posMax as lower and upper bounds
// Done By: Brandon Ting En Junn 2101751
//------------------------------------------------------------------------------------------------------------------------------
void UpdatePosition()
{
	for (int i = 0; i < pSize; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			double position = particlePosition[i][j] + particleVelocity[i][j];

			// Position Limits
			if (position < posMin) {
				particlePosition[i][j] = posMin;
			} else if (position > posMax) {
				particlePosition[i][j] = posMax;
			} else {
				particlePosition[i][j] = position;
			}
		}
	}
}
//------------------------------------------------------------------------------------------------------------------------------

int main()
{
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	srand(time(0));

	// Directory e.g. PSO
	string psoDirectory = "PSO";
	CreateDirectory(psoDirectory.c_str(), NULL);

	string psoInfo = psoDirectory + "\\\\" + psoDirectory + ".txt";
	ofstream outfileo1Info(psoInfo.c_str(), ios::trunc);
	cout << "PSO with constant w and vc" << endl;
	outfileo1Info << "PSO with constant w and vc" << endl << endl;
	outfileo1Info << "Inertia Weight (w): " << w << endl;
	outfileo1Info << "Velocity Clamping (vc): Yes" << endl;
	outfileo1Info << "Constant c (c1): " << c1 << endl;
	outfileo1Info << "Constant c (c2): " << c2 << endl;
	outfileo1Info.close();
	
	//---------------------------------------------------------------------------------------------------------------------------
	// Benchmark Automation
	// Done By: Brandon Ting En Junn 2101751
	//---------------------------------------------------------------------------------------------------------------------------
	for (; BENCHMARK < 11; BENCHMARK++)
	{
		initialiseBenchmark_Range();

		// Directory e.g. PSO/PSO-Sphere
		string benchmarkDirectory = psoDirectory + "\\\\" + psoDirectory + benchmarkFunction;
		CreateDirectory(benchmarkDirectory.c_str(), NULL);
		
		//---------------------------------------------------------------------------------------------------------------------------
		// Experiment Automation
		// Done By: Brandon Ting En Junn 2101751
		//---------------------------------------------------------------------------------------------------------------------------
		for (int i = 0; i < EXPERIMENT; i++)
		{
			// File e.g. PSO/PSO-Sphere/PSO-Sphere-Result1.txt
			string outfile1 = benchmarkDirectory + "\\\\" + psoDirectory + benchmarkFunction + "-Result" + to_string(i + 1) + ".txt";
			ofstream outfileo1(outfile1.c_str(),ios::trunc);
			
			cout << "Experiment " << i + 1 << "..."<< endl;
		
			//CPU Time
			clock_t start, end;
			start = clock();
	//		srand(time(0));
	    
			//---------------------------------------------------------------------------------------------------------------------------
			// Generate Population
			// Done By: Brandon Ting En Junn 2101751
			//---------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
			SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
			cout << "Generate Population..." << endl;
			SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif

	   		for(int i = 0 ; i < pSize ; i++)
	   		{
	      		for(int j = 0 ; j < dimension ; j++)
	      		{
	         		r = getrandom(-rangeMin,rangeMax);
	         		r = r / rangeDiv;
	         		particlePosition[i][j] = r;
	         		particlePBest[i][j] = r;
	      		}

#if DEMO == 1
		      	if (i == 0 || i == pSize - 1) {
		      		cout<<"Particle "<<i+1<<" Position"<<endl;
			    	for(int j = 0 ; j < dimension ; j++)
			    	{
			        	cout<<fixed<<setprecision(3)<<particlePosition[i][j]<<"\t";
			      	}      
			      	cout<<endl<<endl;
				}
#endif
	   		}

#if DEMO == 1
			for (int i = 0; i < pSize; i++)
			{
				if (i == 0 || i == pSize - 1) {
				    cout<<"Particle "<<i+1<<" Velocity"<<endl;
					for(int j = 0 ; j < dimension ; j++)
					{
					    cout<<fixed<<setprecision(3)<<particleVelocity[i][j]<<"\t";
					}      
					cout<<endl<<endl;	
				}	
			} 		
#endif
			// getch();
			//---------------------------------------------------------------------------------------------------------------------------

			//---------------------------------------------------------------------------------------------------------------------------
			// Fitness Evaluation 1
			// Done By: Brandon Ting En Junn 2101751
			//---------------------------------------------------------------------------------------------------------------------------
			for(int i = 0 ; i < pSize ; i++)
			{
				fit[i] = Fitness(particlePosition[i]);
			    // cout<<fixed<<setprecision(3)<<fit[i]<<endl;
				sumFit = 0;
			} 

			//---------------------------------------------------------------------------------------------------------------------------
			// Determine PBest and GBest
			// Done By: Brandon Ting En Junn 2101751
			//---------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
			SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
			cout << "Initial PBest and GBest..." << endl;
			SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif
			//PBest
			for(int i = 0 ; i < pSize ; i++)
			{
				particlePBestFV[i] = fit[i];

#if DEMO == 1
		      	if (i == 0 || i == pSize - 1) {
		      		cout<<"Particle "<<i+1<<" PBest: "<<particlePBestFV[i]<<endl;
			    	for(int j = 0 ; j < dimension ; j++)
			    	{
			        	cout<<fixed<<setprecision(3)<<particlePBest[i][j]<<"\t";
			      	}      
			      	cout<<endl<<endl;
				}
#endif
			}
			
			//GBest
			particleGBestFV = numeric_limits<double>::max();
			for (int i = 0; i < pSize; i++)
			{
				if (particlePBestFV[i] < particleGBestFV)
				{
					particleGBestFV = particlePBestFV[i];
					GBestIndex = i;
				}
			}

#if DEMO == 1
			cout<<"Particle "<<GBestIndex+1<<" GBest: "<<particleGBestFV<<endl;
#endif
			for (int j = 0; j < dimension; j++)
			{
				particleGBest[j] = particlePBest[GBestIndex][j];
#if DEMO == 1
				cout<<fixed<<setprecision(3)<<particleGBest[j]<<"\t";
#endif
			}

#if DEMO == 1
			cout << endl << endl;
#endif

			// cout<<endl;
			// outfileo1<<"Gen\tMinimum"<<endl;
			//---------------------------------------------------------------------------------------------------------------------------
			
			//---------------------------------------------------------------------------------------------------------------------------
			// Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)]
			// Done By: Brandon Ting En Junn 2101751
			//---------------------------------------------------------------------------------------------------------------------------
			for(int i = 0 ; i < gen ; i++)
			{
	      		//------------------------------------------------------------------------------------------------------------------------------
				// Update Velocity
				// Done By: Brandon Ting En Junn 2101751
				//------------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
				SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
				cout << "Update Velocity..." << endl;
				SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif

				UpdateVelocity();

#if DEMO == 1
				for (int i = 0; i < pSize; i++)
				{
					if (i == 0 || i == pSize - 1) {
				      	cout<<"Particle "<<i+1<<" Velocity"<<endl;
					    for(int j = 0 ; j < dimension ; j++)
					    {
					        cout<<fixed<<setprecision(3)<<particleVelocity[i][j]<<"\t";
					    }      
					    cout<<endl<<endl;	
					}	
				}
#endif
				//------------------------------------------------------------------------------------------------------------------------------

				//------------------------------------------------------------------------------------------------------------------------------
				// Update Position
				// Done By: Brandon Ting En Junn 2101751
				//------------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
				SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
				cout << "Update Position..." << endl;
				SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif

				UpdatePosition();

#if DEMO == 1
				for (int i = 0; i < pSize; i++)
				{
					if (i == 0 || i == pSize - 1) {
				      	cout<<"Particle "<<i+1<<" Position"<<endl;
					    for(int j = 0 ; j < dimension ; j++)
					    {
					        cout<<fixed<<setprecision(3)<<particlePosition[i][j]<<"\t";
					    }      
					    cout<<endl<<endl;	
					}	
				}
#endif
				//------------------------------------------------------------------------------------------------------------------------------

				//------------------------------------------------------------------------------------------------------------------------
				// Fitness Evaluation 2
				// Done By: Brandon Ting En Junn 2101751
				//------------------------------------------------------------------------------------------------------------------------
				for(int i = 0 ; i < pSize ; i++)
				{
					fit[i] = Fitness(particlePosition[i]);
			      	// cout<<fixed<<setprecision(3)<<fit[i]<<endl;
					sumFit = 0;
				} 
				//------------------------------------------------------------------------------------------------------------------------

	      		//------------------------------------------------------------------------------------------------------------------------
				// Update PBest and GBest
				// Done By: Brandon Ting En Junn 2101751
				//------------------------------------------------------------------------------------------------------------------------
#if DEMO == 1
				SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);
				cout << "Updated PBest and GBest..." << endl;
				SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#endif
				//PBest
				for(int i = 0 ; i < pSize ; i++)
				{
					if(fit[i] < particlePBestFV[i])
					{
						particlePBestFV[i] = fit[i];

						for (int j = 0; j < dimension; j++)
						{
							particlePBest[i][j] = particlePosition[i][j];
						}
					}
				}

#if DEMO == 1
				for (int i = 0; i < pSize; i++)
				{
					if (i == 0 || i == pSize - 1) {
				      	cout<<"Particle PBest "<<i+1<<": "<<particlePBestFV[i]<<endl;
					    for(int j = 0 ; j < dimension ; j++)
					    {
					        cout<<fixed<<setprecision(3)<<particlePBest[i][j]<<"\t";
					    }      
					    cout<<endl<<endl;	
					}	
				}
#endif

				//GBest
				for (int i = 0; i < pSize; i++)
				{
					if (particlePBestFV[i] < particleGBestFV)
					{
						particleGBestFV = particlePBestFV[i];
						GBestIndex = i;
					}
				}
				// cout<<"The group best particle is personal best particle "<<GBestIndex+1<<" with fitness of "<<particleGBestFV<<endl;
				for (int j = 0; j < dimension; j++)
				{
					particleGBest[j] = particlePBest[GBestIndex][j];
					// cout<<particleGBest[j]<<" ";
				}
				// cout<<endl;
		
		      	// fit1 = 0;
		      	// fit2 = 0;

#if DEMO == 1
				cout<<"Particle "<<GBestIndex<<" GBest: "<<particleGBestFV<<endl;
				for(int j = 0 ; j < dimension ; j++)
				{
					cout<<fixed<<setprecision(3)<<particleGBest[j]<<"\t";
				}      
				cout<<endl<<endl;
#endif
		       
		      	outfileo1<<setprecision(6)<<particleGBestFV<<endl;
		      	// cout<<fixed<<setprecision(3)<<particleGBestFV<<endl;
		      	// getch();

#if DEMO == 1
				cout << "SUCCESSFULLY COMPLETED DEMO...\nPress Any Key to Exit..." << endl;
				getch();
				exit(0);
#endif
	    	} //Termination Criteria (Maximum Generation) [Can modify starting from here (GA, DE, PSO)] LOOP
	    	//------------------------------------------------------------------------------------------------------------------------ 

	    	outfileo1<<endl<<endl;
	    	cout<<"Result"<<endl;
	    	
	    	//Output Final Generation Lowest
	//    	cout<<fixed<<setprecision(3)<<lFv<<" "<<lFvIndex+1<<endl<<endl;
	//	    for(int j = 0 ; j < dimension ; j++)
	//	    {
	//	       cout<<fixed<<setprecision(3)<<chromosome[lFvIndex][j]<<"\t";
	//	       outfileo1<<setprecision(6)<<chromosome[lFvIndex][j]<<"\n";
	//	    }
			
			//Output Group Best in Experiment
	    	cout<<fixed<<setprecision(3)<<particleGBestFV<<endl<<endl;
	    	for(int j = 0 ; j < dimension ; j++)
	    	{
	    		cout<<fixed<<setprecision(3)<<particleGBest[j]<<"\t";
	       		outfileo1<<setprecision(6)<<particleGBest[j]<<"\n";
	    	}
		    
		    cout<<endl;
		    outfileo1<<endl;
		    end = clock();
		    cout<<"Time required for execution: "<< (double)(end-start)/CLOCKS_PER_SEC<<" seconds."<<"\n\n";  
		    outfileo1<<(double)(end-start)/CLOCKS_PER_SEC<<"\n\n";
		    outfileo1.close();
			cout<<endl<<endl;
		 
			resetExperiment();
		} //Experiment Automation LOOP
		//------------------------------------------------------------------------------------------------------------------------

	} //Benchmark Automation LOOP
	//------------------------------------------------------------------------------------------------------------------------

	cout << "SUCCESSFULLY COMPLETED...\nPress Any Key to Exit..." << endl;
	getch();
	
	return 0;   
}