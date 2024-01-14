#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>

typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;

} Point;

typedef struct {
	double dist;
	int index;
} MinDist;

//Reads the image dimensions a x b pixels
void readImageSize(FILE* ifp, int* K, int* a, int* b)
{
	fscanf(ifp, "%d\n", K);
	printf("%d\n", *K);

	fscanf(ifp, "%d\n", a);
	printf("%d\n", *a);

	fscanf(ifp, "%d\n", b);
	printf("%d\n", *b);
}

//reads the ifp file and stores in structure
void readPoints(FILE* ifp, Point* points, int num_points)
{
	int i;
	for (i = 0;i < num_points;i++)
	{
		fscanf(ifp, "%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
	}
}

////Initialize random points as assumed means
//void initialize(Point* mean, int K, int num_points, Point* points)
//{
//	int i, a, p = 2;
//	srand(time(NULL));
//	for (i = 0;i < K;i++)
//	{
//		a = num_points / p;
//
//		mean[i]._r = points[a]._r;
//		mean[i]._g = points[a]._g;
//		mean[i]._b = points[a]._b;
//		mean[i]._m = points[a]._m;
//		mean[i]._n = points[a]._n;
//		p++;
//	}
//}

void initialize(Point* mean, int K, int num_points, Point* points)
{
	int i, a, p = 2;
	unsigned int randVal;

	#pragma omp parallel private(i, a, randVal)
		{
	#pragma omp for
			for (i = 0; i < K; i++)
			{
				// Generisanje slučajnog broja uz pomoć rand_s
				rand_s(&randVal);
				a = randVal % num_points;

				mean[i]._r = points[a]._r;
				mean[i]._g = points[a]._g;
				mean[i]._b = points[a]._b;
				mean[i]._m = points[a]._m;
				mean[i]._n = points[a]._n;
			}
		}
}

//All points having no clusters
int IntClusterMem(int* cluster, int num_points)
{
	int i;
	for (i = 0;i < num_points;i++)
	{
		cluster[i] = -1;
	}
}

//Distance
double calculateDistance(Point point1, Point point2)
{
	return sqrt((pow((point1._r - point2._r), 2) + pow((point1._g - point2._g), 2) + pow((point1._b - point2._b), 2) + pow((point1._m - point2._m), 2) + pow((point1._n - point2._n), 2)));
}

//double calculateDistance(Point point1, Point point2) {
//	double dr = point1._r - point2._r;
//	double dg = point1._g - point2._g;
//	double db = point1._b - point2._b;
//	double dm = point1._m - point2._m;
//	double dn = point1._n - point2._n;
//
//	return sqrt(dr * dr + dg * dg + db * db + dm * dm + dn * dn);
//}


//to calculate which cluster is the point belonging to.
//int pointsCluster(Point point, Point* mean, int K)
//{
//	int parent = 0;
//	double dist = 0;
//	double minDist = calculateDistance(point, mean[0]);
//	int i;
//	for (i = 1;i < K;i++)
//	{
//		dist = calculateDistance(point, mean[i]);
//		if (minDist >= dist)
//		{
//			parent = i;
//			minDist = dist;
//		}
//	}
//	return parent;
//}

int pointsCluster(Point point, Point* mean, int K) {
	MinDist minDist = { calculateDistance(point, mean[0]), 0 };
	int i;

	#pragma omp parallel private(i)
		{
			MinDist localMinDist = minDist;
			#pragma omp for nowait
				for (i = 1; i < K; i++) {
					double dist = calculateDistance(point, mean[i]);
					if (dist < localMinDist.dist) {
						localMinDist.dist = dist;
						localMinDist.index = i;
					}
				}
			#pragma omp critical
				{
					if (localMinDist.dist < minDist.dist) {
						minDist = localMinDist;
					}
				}
		}

	return minDist.index;
}


//calculate new mean
//void calcNewMean(Point* points, int* cluster, Point* mean, int K, int num_points)
//{
//	Point* newMean = malloc(sizeof(Point) * K);
//	int* members = malloc(sizeof(int) * K);
//	int i;
//	for (i = 0;i < K;i++)
//	{
//		members[i] = 0;
//		newMean[i]._r = 0;
//		newMean[i]._g = 0;
//		newMean[i]._b = 0;
//		newMean[i]._m = 0;
//		newMean[i]._n = 0;
//	}
//	for (i = 0;i < num_points;i++)
//	{
//		members[cluster[i]]++;
//		newMean[cluster[i]]._r += points[i]._r;
//		newMean[cluster[i]]._g += points[i]._g;
//		newMean[cluster[i]]._b += points[i]._b;
//		newMean[cluster[i]]._m += points[i]._m;
//		newMean[cluster[i]]._n += points[i]._n;
//	}
//	for (i = 0;i < K;i++)
//	{
//		if (members[i] != 0.0)
//		{
//			newMean[i]._r /= members[i];
//			newMean[i]._g /= members[i];
//			newMean[i]._b /= members[i];
//			newMean[i]._m /= members[i];
//			newMean[i]._n /= members[i];
//		}
//		else
//		{
//			newMean[i]._r = 0;
//			newMean[i]._g = 0;
//			newMean[i]._b = 0;
//			newMean[i]._m = 0;
//			newMean[i]._n = 0;
//		}
//	}
//	for (i = 0;i < K;i++)
//	{
//		mean[i]._r = newMean[i]._r;
//		mean[i]._g = newMean[i]._g;
//		mean[i]._b = newMean[i]._b;
//		mean[i]._m = newMean[i]._m;
//		mean[i]._n = newMean[i]._n;
//	}
//}

void calcNewMean(Point* points, int* cluster, Point* mean, int K, int num_points) {
	Point* newMeanTotal = (Point*)calloc(K, sizeof(Point));
	int* membersTotal = (int*)calloc(K, sizeof(int));
	int i;

	// Paralelizacija akumulacije sume i brojanja članova
	#pragma omp parallel private(i)
		{
			Point* localMean = (Point*)calloc(K, sizeof(Point));
			int* localMembers = (int*)calloc(K, sizeof(int));

			// Paralelno sumiranje
		#pragma omp for
			for (i = 0; i < num_points; i++) {
				int idx = cluster[i];
				localMean[idx]._r += points[i]._r;
				localMean[idx]._g += points[i]._g;
				localMean[idx]._b += points[i]._b;
				localMean[idx]._m += points[i]._m;
				localMean[idx]._n += points[i]._n;
				localMembers[idx]++;
			}

			// Kombinovanje rezultata iz svih niti
		#pragma omp critical
			{
				for (int i = 0; i < K; i++) {
					newMeanTotal[i]._r += localMean[i]._r;
					newMeanTotal[i]._g += localMean[i]._g;
					newMeanTotal[i]._b += localMean[i]._b;
					newMeanTotal[i]._m += localMean[i]._m;
					newMeanTotal[i]._n += localMean[i]._n;
					membersTotal[i] += localMembers[i];
				}
			}

			free(localMean);
			free(localMembers);
		}

		// Izračunavanje novih srednjih vrednosti
		for (int i = 0; i < K; i++) {
			if (membersTotal[i] > 0) {
				mean[i]._r = newMeanTotal[i]._r / membersTotal[i];
				mean[i]._g = newMeanTotal[i]._g / membersTotal[i];
				mean[i]._b = newMeanTotal[i]._b / membersTotal[i];
				mean[i]._m = newMeanTotal[i]._m / membersTotal[i];
				mean[i]._n = newMeanTotal[i]._n / membersTotal[i];
			}
		}

	free(newMeanTotal);
	free(membersTotal);
}

//int chkConvrg(int* before_clusters, int* after_cluster, int num_points, float tol)
//{
//	int i;
//	tol = num_points * tol;
//	for (i = 0;i < num_points;i++)
//		if (abs(before_clusters[i] - after_cluster[i]) > tol)
//			return -1;
//	return 0;
//}

int chkConvrg(int* before_clusters, int* after_cluster, int num_points, float tol) {
	int i;
	int result = 0; // 0 za konvergenciju, 1 za nekonvergenciju
	tol = num_points * tol;

	#pragma omp parallel for reduction(|:result)
		for (i = 0; i < num_points; i++) {
			if (abs(before_clusters[i] - after_cluster[i]) > tol) {
				result = 1; // Ako postoji nekonvergencija, postavi na 1
				// Ne možemo odmah prekinuti petlju, ali zabeležimo nekonvergenciju
			}
		}

	return (result == 0) ? 0 : -1;
}


int main(int argc, char* argv[])
{
	int iter = 0;
	int K;
	int num_points;
	int i;
	int job_done = 0;
	int x, y;
	float tol;

	LARGE_INTEGER frequency;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	double interval;


	Point* mean;
	Point* points;
	Point* get_points;
	int* formed_clusters;
	int* before_clusters;
	int* after_cluster;
	int thread_no;

	printf("Enter tolerence: ");
	scanf("%f", &tol);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);

	thread_no = atoi(argv[3]);

	FILE* ifp;
	ifp = fopen(argv[1], "r");
	readImageSize(ifp, &K, &x, &y);
	num_points = x * y;
	points = (Point*)malloc(sizeof(Point) * num_points);
	readPoints(ifp, points, num_points);
	fclose(ifp);

	before_clusters = (int*)malloc(sizeof(int) * num_points);
	after_cluster = (int*)malloc(sizeof(int) * num_points);
	mean = malloc(sizeof(Point) * K);

	//initializing to default values
	initialize(mean, K, num_points, points);
	IntClusterMem(before_clusters, num_points);
	IntClusterMem(after_cluster, num_points);


	while (1)
	{
		iter++;
		printf("Iteration %d\n", iter);

		omp_set_num_threads(thread_no);

	#pragma omp parallel for default(shared) 

		for (i = 0;i < num_points;i++)
		{

			after_cluster[i] = pointsCluster(points[i], mean, K);
		}

		calcNewMean(points, after_cluster, mean, K, num_points);
		printf("New Centroids are calculated!\n");

		if (chkConvrg(after_cluster, before_clusters, num_points, tol) == 0)
		{

			printf("K-mean algorithm Converged!\n");
			job_done = 1;
		}
		else
		{
			printf("Not converged!\n");
			for (i = 0;i < num_points;i++)
				before_clusters[i] = after_cluster[i];
		}

		if (job_done == 1)
			break;

	}


	FILE* ofp = fopen(argv[2], "w");

	for (i = 0;i < K;i++)
		fprintf(ofp, "%d,%d,%d,%d,%d\n", (int)mean[i]._r, (int)mean[i]._g, (int)mean[i]._b, (int)mean[i]._m, (int)mean[i]._n);

	for (i = 0;i < num_points;i++)
		fprintf(ofp, "%d,%d,%d,%d,%d,%d\n", (int)points[i]._r, (int)points[i]._g, (int)points[i]._b, (int)points[i]._m, (int)points[i]._n, after_cluster[i] + 1);

	fclose(ofp);

	QueryPerformanceCounter(&end);
	interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

	printf("Total time elapsed in forming clusters : %.8f sec\n", interval);

	return 0;
}
