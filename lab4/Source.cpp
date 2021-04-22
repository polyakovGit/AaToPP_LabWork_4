#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <ctime>

void Task1ParallelFor(double** arr, long long N, int k, double& avg) {
#pragma omp parallel num_threads(k)
	{
#pragma omp for
		for (long long i = 0; i < N; ++i) {
			for (long long j = 0; j < N; ++j) {
				arr[i][j] = sin(i) + cos(j);
				avg += arr[i][j];
			}
		}
	}
avg /= N*N;
}

int main() {
	const long long N = 20;
	double avg = 0;
	clock_t start, end;
	printf("%lld\n", N);
	double** arr = new double* [N];
	for (int i = 0; i < N; ++i)
		arr[i] = new double[N];
	start = clock();
	for (long long i = 0; i < N; ++i)
		for (long long j = 0; j < N; ++j) {
			arr[i][j] = sin(i) + cos(j);
			avg += arr[i][j];
		}
	avg /= N*N;
	end = clock();
	printf("SERIAL value: %f, time:%f\n", avg, (double)end - start / CLOCKS_PER_SEC);

	const int arrThreadSize = 3;
	int arrThreads[arrThreadSize] = { 2,4,8 };
	for (int k = 0; k < arrThreadSize; ++k)
	{
		avg = 0;
		start = clock();
		Task1ParallelFor(arr, N, arrThreads[k], avg);
		end = clock();
		printf("Task1ParallelFor value: %f, time: %f, threads: %d\n", avg, (double)end - start / CLOCKS_PER_SEC, arrThreads[k]);
	}
	return 0;
}