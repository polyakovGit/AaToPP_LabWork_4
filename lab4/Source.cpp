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
	const long long N = 10000;
	double serialAvg = 0;
	clock_t start, end;
	printf("%lld\n", N);
	double** arr = new double* [N];
	for (int i = 0; i < N; ++i)
		arr[i] = new double[N];
	start = clock();
	for (long long i = 0; i < N; ++i)
		for (long long j = 0; j < N; ++j) {
			arr[i][j] = sin(i) + cos(j);
			serialAvg += arr[i][j];
		}
	serialAvg /= N * N;
	end = clock();
	printf("SERIAL value: %f, time:%f\n", serialAvg, (double)(end - start) / CLOCKS_PER_SEC);

	Func arrFunc[4] = { {Task1ParallelFor,"Task1ParallelFor"},{Task1Reduction,"Task1Reduction"},
		{Task1Atomic,"Task1Atomic"}, {Task1Critical,"Task1Critical"} };
	//void((*ptrFunc[4]))(double**, long long, int, double&);
	//ptrFunc[0] = Task1ParallelFor;
	//ptrFunc[1] = Task1Reduction;
	//ptrFunc[2] = Task1Atomic;
	//ptrFunc[3] = Task1Critical;
	const int arrThreadSize = 3;
	int arrThreads[arrThreadSize] = { 2,4,8 };
	for (int k = 0; k < arrThreadSize; ++k)
	{
		for (int i = 0; i < 4; ++i) {
			double time = 0;
			for (int j = 0; j < 1; ++j) {
				double parallelAvg=0;
				start = clock();
				arrFunc[i].ptrFunc(arr, N, arrThreads[k], parallelAvg);
				//ptrFunc[i](arr, N, arrThreads[k], parallelAvg);
				end = clock();
				time += (double)(end - start) / CLOCKS_PER_SEC;
				printf("%f ", parallelAvg);
			}
			time /= 1;
			printf("Func: %s, average time:%f threads: %d\n", arrFunc[i].name.c_str(),time, arrThreads[k]);
		}


		//double parallelAvg = 0;
		//start = clock();
		//Task1ParallelFor(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1ParallelFor value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
		//parallelAvg = 0;
		//start = clock();
		//Task1Reduction(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1Reduction value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
		//parallelAvg = 0;
		//start = clock();
		//Task1Atomic(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1Atomic value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
		//parallelAvg = 0;
		//start = clock();
		//Task1Critical(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1Critical value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
	}
	return 0;
}