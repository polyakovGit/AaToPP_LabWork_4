#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>

struct Func {
	void((*ptrFunc))(double**, long long, int, double&);
	std::string name;
};

void Task1ParallelFor(double** arr, long long N, int k, double& avg) {
#pragma omp parallel for num_threads(k) 
	for (long long i = 0; i < N; ++i) {
		double sum = 0;
		for (long long j = 0; j < N; ++j) {
			arr[i][j] = sin(i) + cos(j);
			sum += arr[i][j];
		}
		avg += sum;
	}
	avg /= (N * N);
}

void Task1Reduction(double** arr, long long N, int k, double& avg) {
	double sum = 0;
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
		long long currThread = omp_get_thread_num();
		for (long long i = currThread * N / k; i < (currThread + 1) * N / k; ++i) {
			for (long long j = 0; j < N; ++j) {
				arr[i][j] = sin(i) + cos(j);
				sum += arr[i][j];
			}
		}
	}
	avg = sum / (N * N);
}


void Task1Atomic(double** arr, long long N, int k, double& avg) {
	double sum = 0;
#pragma omp parallel num_threads(k)
	{
		long long currThread = omp_get_thread_num();
		for (long long i = currThread * N / k; i < (currThread + 1) * N / k; ++i) {
			for (long long j = 0; j < N; ++j) {
				arr[i][j] = sin(i) + cos(j);
#pragma omp atomic
				sum += arr[i][j];
			}
		}
	}
	avg = sum / (N * N);
}

void Task1Critical(double** arr, long long N, int k, double& avg) {
	double sum = 0;
#pragma omp parallel num_threads(k)
	{
		long long currThread = omp_get_thread_num();
		for (long long i = currThread * N / k; i < (currThread + 1) * N / k; ++i) {
			for (long long j = 0; j < N; ++j) {
				arr[i][j] = sin(i) + cos(j);
#pragma omp critical
				sum += arr[i][j];
			}
		}
	}
	avg = sum / (N * N);
}

int main() {
	//const long long N = 10000;
	//double serialAvg = 0;
	//clock_t start, end;
	//printf("%lld\n", N);
	//double** arr = new double* [N];
	//for (int i = 0; i < N; ++i)
	//	arr[i] = new double[N];
	//start = clock();
	//for (long long i = 0; i < N; ++i)
	//	for (long long j = 0; j < N; ++j) {
	//		arr[i][j] = sin(i) + cos(j);
	//		serialAvg += arr[i][j];
	//	}
	//serialAvg /= N * N;
	//end = clock();
	//printf("SERIAL value: %f, time:%f\n", serialAvg, (double)(end - start) / CLOCKS_PER_SEC);

	//Func arrFunc[4] = { {Task1ParallelFor,"Task1ParallelFor"},{Task1Reduction,"Task1Reduction"},
	//	{Task1Atomic,"Task1Atomic"}, {Task1Critical,"Task1Critical"} };
	//const int arrThreadSize = 3;
	//int arrThreads[arrThreadSize] = { 2,4,8 };
	//for (int k = 0; k < arrThreadSize; ++k)
	//{
	//	for (int i = 0; i < 1; ++i) {
	//		double time = 0;
	//		for (int j = 0; j < 1; ++j) {
	//			double parallelAvg = 0;
	//			start = clock();
	//			arrFunc[i].ptrFunc(arr, N, arrThreads[k], parallelAvg);
	//			//ptrFunc[i](arr, N, arrThreads[k], parallelAvg);
	//			end = clock();
	//			time += (double)(end - start) / CLOCKS_PER_SEC;
	//			printf("%f ", parallelAvg);
	//		}
	//		time /= 1;
	//		printf("Func: %s, average time:%f threads: %d\n", arrFunc[i].name.c_str(), time, arrThreads[k]);
	//	}
	//}

	const int SIZE = 5;
	int** arrMulA = new int* [SIZE];
	int** arrMulB = new int* [SIZE];
	int** arrMulC = new int* [SIZE];
	for (int i = 0; i < SIZE; ++i) {
		arrMulA[i] = new int[SIZE];
		arrMulB[i] = new int[SIZE];
		arrMulC[i] = new int[SIZE];
	}

	for (int i = 0; i < SIZE; ++i)
		for (int j = 0; j < SIZE; ++j)
		{
			arrMulA[i][j] = arrMulB[i][j] = 2;
		}

#pragma omp parallel for
	for (int i = 0; i < SIZE; ++i) {
		for (int j = 0; j < SIZE; ++j) {
			int temp = 0;
			for (int k = 0; k < SIZE; ++k)
				temp += arrMulA[i][k] * arrMulB[k][j];
			arrMulC[i][j] = temp;
		}
	}
	for (int i = 0; i < SIZE; ++i) {
		for (int j = 0; j < SIZE; ++j)
			printf("%d ", arrMulC[i][j]);
		printf("\n");
	}

	return 0;
}