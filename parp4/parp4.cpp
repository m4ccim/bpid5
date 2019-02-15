#include "pch.h"
#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdio.h>
#include <intrin.h>
#include <WinSock2.h>

using namespace std;




double getMadhava(int iterations) {
	double pi = 1;
	for (int i = 1; i < iterations; i++) {
		if (i % 2 == 0)
			pi += (double)1 / ((i * 2 + 1)*pow(3, i));
		else if (i % 2 == 1)
			pi -= (double)1 / ((i * 2 + 1)*pow(3, i));
	}
	return pi * sqrt(12);
}

double getLeibnizSum(long long iterations) {
	double pi = 1;
	for (long long i = 1; i < iterations; i++) {
		if (i % 2 ==  0)
			pi += (double)1 / (i * 2 + 1);
		else if (i % 2 == 1)
			pi -= (double)1 / (i * 2 + 1);
	}	return pi * 4;
}

double getBBP(int digit) {
	double pi = 0.0;
	for (int k =0; k<digit; ++k)
	{
		pi += ((double)(double)4 / (8 * k + 1) 
			- (double)2 / (8 * k + 4) 
			- (double)1 / (8 * k + 5) 
			- (double)1 / (8 * k + 6)) 
			/ pow(16, k);
	}
	return pi;
}

double getBellard(int iterations) {
	double pi = 0.0;
	for (int n = 0; n < iterations; n++) {
		pi += pow((-1), n) / pow(2, 10 * n)
			*(-pow(2, 5) / (4 * n + 1) - (double)1 / (4 * n + 3)
				+ pow(2, 8) / (10 * n + 1) - pow(2, 6) / (10 * n + 3)
				- pow(2, 2) / (10 * n + 5) - pow(2, 2) / (10 * n + 7) 
				+ (double)1 / (10 * n + 9));
	}
	return pi / pow(2, 6);
}

double getBellardParallel(int iterations) {
	double pi = 0.0;

#pragma omp parallel for reduction (+:pi)
	for (int n = 0; n < iterations; n++) {
		double P[6];
		double Z[7];

		P[0] = pow(-1, n);
		P[1] = pow(2, 10 * n);
		P[2] = pow(2, 5); 
		P[3] = pow(2, 8);
		P[4] = pow(2, 6); 
		P[5] = pow(2, 2); 

		Z[0] = 4 * n + 1; 
		Z[1] = 4 * n + 3; 
		Z[2] = 10 * n + 1; 
		Z[3] = 10 * n + 3; 
		Z[4] = 10 * n + 5; 
		Z[5] = 10 * n + 7; 
		Z[6] = 10 * n + 9; 

		pi += P[0] / P[1] * (-P[2] / Z[0] - (double)1 
			/ Z[1]+ P[3] / Z[2] - P[4] / Z[3] - P[5]
			/ Z[4] - P[5] / Z[5]+ (double)1 / Z[6]);
	}
	return pi / pow(2, 6);
}

typedef struct  {
	double piPart;
	int begin;
	int end;
} iteration;

DWORD WINAPI func(LPVOID a) {
	double P[6];
	double Z[7];
	iteration * itr = (iteration *)a;
	double pi = 0;
	for (int n = itr->begin; n < itr->end; n++) {

		P[0] = pow(-1, n);
		P[1] = pow(2, 10 * n);
		P[2] = pow(2, 5);
		P[3] = pow(2, 8);
		P[4] = pow(2, 6);
		P[5] = pow(2, 2);

		Z[0] = 4 * n + 1;
		Z[1] = 4 * n + 3;
		Z[2] = 10 * n + 1;
		Z[3] = 10 * n + 3;
		Z[4] = 10 * n + 5;
		Z[5] = 10 * n + 7;
		Z[6] = 10 * n + 9;

		pi += P[0] / P[1] * (-P[2] / Z[0] - (double)1 
			/ Z[1] + P[3] / Z[2] - P[4] / Z[3] - P[5] 
			/ Z[4] - P[5] / Z[5] + (double)1 / Z[6]);
	}
	itr->piPart = pi;
	return 1;
}

double getBellardParallelWin(int iterations) {
	int Regs[4];
	__cpuid(Regs, 0x00000001);
	int CoresCount = (Regs[1] >> 16) & 0xFF;
	double pi = 0.0;
	
	iteration* itr = new iteration[CoresCount];

	for (int n = 0; n < CoresCount; n++) {
		itr[n].begin = n* (iterations / CoresCount);
		itr[n].end = (n + 1)* (iterations / CoresCount);
	}

	HANDLE* thrd = new HANDLE[CoresCount];

	for (int n = 0; n < CoresCount;n++) {

			thrd[n] = CreateThread(0, 0, func, &itr[n], 0, 0);
	}

	WaitForMultipleObjects(CoresCount, thrd, true, INFINITE);

	for (int i = 0; i < CoresCount; i++) {
		pi += itr[i].piPart;
		CloseHandle(thrd[i]);
	}

	return pi / pow(2, 6);
}




double round(double x, int n) {
	return floor(x * pow(10, n) + .5) / pow(10,n);
}

int main()
{

#ifdef _OPENMP
	printf("_OPENMP Defined\n");
#else
	printf("_OPENMP UnDefined\n");
#endif
	
	cout << "PI = 3.14159265358979\n";
	cout.precision(15);
	double pi = 0.0;
	for (int i = 0; i < 10; i++) {
		double  dif = 0.0;
		double start = omp_get_wtime(); //start the timer

		pi = getBellardParallel(16000000 +i);
		cout << pi;

		double end = omp_get_wtime();// end the timer
		dif = end - start; // stores the difference in dif
		printf(" Spent time: %lf seconds\n", dif);
	}
	cout << "My PI= " << pi;
	
}
