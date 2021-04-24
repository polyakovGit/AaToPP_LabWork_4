#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <string>
#include <fstream>

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
				{
					sum += arr[i][j];
				}
			}
		}
	}
	avg = sum / (N * N);
}

void ParallelMul(int** arrMulA, int** arrMulB, int** arrMulC, int n, int numThreads) {
#pragma omp parallel for num_threads(numThreads) schedule(dynamic,1)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			int temp = 0;
			for (int k = 0; k < n; ++k)
				temp += arrMulA[i][k] * arrMulB[k][j];
			arrMulC[i][j] = temp;
		}
	}
}

void FillArrays(int** arrMulA, int** arrMulB, int n) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			arrMulA[i][j] = arrMulB[i][j] = j;//без рандома
}

void quickSort(double arr[], int n) {
	int i = 0, j = n - 1;
	double pivot = arr[n / 2];
	do {
		while (arr[i] < pivot)++i;
		while (arr[j] > pivot)--j;
		if (i <= j) {
			std::swap(arr[i], arr[j]);
			++i;
			--j;
		}
	} while (i <= j);
	if (n > i)quickSort(arr + i, n - i);
	if (j > 0)quickSort(arr, j + 1);
}

void GetAvgTime(double arrTime[], int n, double& avgTime) {
	for (int s = 5; s < 15; ++s)
		avgTime += arrTime[s];
	avgTime /= 10;
}

int main() {
	const long long N = 1000;
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

	const int arrFuncSize = 4;
	Func arrFunc[arrFuncSize] = { {Task1ParallelFor,"Task1ParallelFor"},{Task1Reduction,"Task1Reduction"},
		{Task1Atomic,"Task1Atomic"}, {Task1Critical,"Task1Critical"} };
	const int Task1arrThreadSize = 3;
	int Task1arrThreads[Task1arrThreadSize] = { 2,4,8 };
	const int timeRepeat = 20;
	for (int k = 0; k < Task1arrThreadSize; ++k)
	{
		for (int i = 0; i < arrFuncSize; ++i) {
			double time = 0;
			for (int j = 0; j < timeRepeat; ++j) {
				double parallelAvg = 0;
				start = clock();
				arrFunc[i].ptrFunc(arr, N, Task1arrThreads[k], parallelAvg);
				end = clock();
				time += (double)(end - start) / CLOCKS_PER_SEC;
				printf("%f\n", parallelAvg);
			}
			time /= timeRepeat;
			printf("Func: %s,\taverage time: %f\tthreads: %d\n", arrFunc[i].name.c_str(), time, Task1arrThreads[k]);
		}
	}

	for (int i = 0; i < N; ++i)
		delete[]arr[i];
	delete[] arr;

	const int arrNumSize = 4;
	int arrNums[arrNumSize] = { 100,1000,5000,10000 };
	const int arrThreadSize = 9;
	int arrThread[arrThreadSize] = { 1,2,4,8,10,16,20,24,30 };
	std::ofstream OutputData;
	OutputData.open("OutputData.csv");

	for (int k = 0; k < arrThreadSize; ++k) {
		printf("\t%d\t", arrThread[k]);
		OutputData << ';' << arrThread[k];
	}

	for (int i = 0; i < arrNumSize; ++i)
	{
		printf("\n%d\t", arrNums[i]);
		OutputData << '\n' << arrNums[i];
		int** arrMulA = new int* [arrNums[i]];
		int** arrMulB = new int* [arrNums[i]];
		int** arrMulC = new int* [arrNums[i]];
		for (int elem = 0; elem < arrNums[i]; ++elem) {
			arrMulA[elem] = new int[arrNums[i]];
			arrMulB[elem] = new int[arrNums[i]];
			arrMulC[elem] = new int[arrNums[i]];
		}

		FillArrays(arrMulA, arrMulB, arrNums[i]);

		for (int k = 0; k < arrThreadSize; ++k)
		{
			double arrTime[timeRepeat];
			//еще цикл для среднего времени 20 итераций взять средние квартили
			for (int t = 0; t < timeRepeat; ++t) {
				start = clock();
				ParallelMul(arrMulA, arrMulB, arrMulC, arrNums[i], arrThread[k]);
				end = clock();
				double time = (double)(end - start) / CLOCKS_PER_SEC;
				arrTime[t] = time;
			}
			quickSort(arrTime, timeRepeat);
			double avgTime = 0;
			GetAvgTime(arrTime, timeRepeat,avgTime);//среднее время по средним квартилям 

			printf("%f\t", avgTime);
			OutputData << ';' << avgTime;
		}
		for (int del = 0; del < arrNums[i]; ++del) {
			delete[] arrMulA[del];
			delete[] arrMulB[del];
			delete[] arrMulC[del];
		}
		delete[] arrMulA;
		delete[] arrMulB;
		delete[] arrMulC;
	}

	return 0;
}