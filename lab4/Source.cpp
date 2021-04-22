#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <ctime>

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
//
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


double arrayFillAndSum(double** arr, long long currI, long long currN, long long N) {
	double sum = 0;
	for (long long i = currI; i < currN; ++i)
		for (long long j = 0; j < N; ++j) {
			arr[i][j] = sin(i) + cos(j);
			sum += arr[i][j];
		}
	return sum;
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
				sum += arr[i][j];//arrayFillAndSum(arr, currThread * N / k, (currThread+1) * N / k, N);
			}
			avg = sum / (N * N);
		}
	}
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

	const int arrThreadSize = 3;
	int arrThreads[arrThreadSize] = { 2,4,8 };
	for (int k = 0; k < arrThreadSize; ++k)
	{
		double parallelAvg = 0;
		//start = clock();
		//Task1ParallelFor(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1ParallelFor value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
		//parallelAvg = 0;
		start = clock();
		Task1Reduction(arr, N, arrThreads[k], parallelAvg);
		end = clock();
		printf("Task1Reduction value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
		//parallelAvg = 0;
		//start = clock();
		//Task1Atomic(arr, N, arrThreads[k], parallelAvg);
		//end = clock();
		//printf("Task1Atomic value: %f, time: %f, threads: %d\n", parallelAvg, (double)(end - start) / CLOCKS_PER_SEC, arrThreads[k]);
	}
	return 0;
}