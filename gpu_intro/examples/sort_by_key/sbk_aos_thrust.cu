#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

// Get the timer value.
static void get_time(double* ret)
{
	volatile struct timespec val;
	clock_gettime(CLOCK_REALTIME, (struct timespec*)&val);
	*ret = (double)0.000000001 * val.tv_nsec + val.tv_sec;
}


void usage(const char* filename)
{
	printf("Sort the random key-value data set of the given length by key.\n");
	printf("Usage: %s <n>\n", filename);
}

using namespace thrust;

static inline pair<float, float> rand_pair()
{
	pair<float, float> val;
	val.first = drand48();
	val.second = drand48();
	return val;
}

int main(int argc, char* argv[])
{
	const int printable_n = 128;

	if (argc != 2)
	{
		usage(argv[0]);
		return 0;
	}

	int n = atoi(argv[1]);
	if (n <= 0)
	{
		usage(argv[0]);
		return 0;
	}

	// Generate the random keys and values on the host.
	host_vector<pair<float, float> > h_pairs(n);
	generate(h_pairs.begin(), h_pairs.end(), rand_pair);

	// Print out the input data if n is small.
	if (n <= printable_n)
	{
		printf("Input data:\n");
		for (int i = 0; i < n; i++)
			printf("(%f, %f)\n", h_pairs[i].first, h_pairs[i].second);
		printf("\n");
	}

	// Strip initialization delay.
	double init_start, init_finish;
	get_time(&init_start);
#ifndef CPU
	int count = 0;
	cudaGetDeviceCount(&count);
#endif
	get_time(&init_finish);

	// Transfer data to the device.
	double load_start, load_finish;
	get_time(&load_start);
#ifndef CPU
	device_vector<pair<float, float> > d_pairs = h_pairs;
#endif
	get_time(&load_finish);

	double sort_start, sort_finish;
	get_time(&sort_start);
#ifdef CPU
	sort(h_pairs.begin(), h_pairs.end());
#else
	sort(d_pairs.begin(), d_pairs.end());
	cudaDeviceSynchronize();
#endif
	get_time(&sort_finish);

	// Transfer data back to host.
	double save_start, save_finish;
	get_time(&save_start);
#ifndef CPU
	copy(d_pairs.begin(), d_pairs.end(), h_pairs.begin());
#endif
	get_time(&save_finish);
	
	printf("Init time = %f sec\n", init_finish - init_start);
	printf("Load time = %f sec\n", load_finish - load_start);
	printf("Sort time = %f sec\n", sort_finish - sort_start);
	printf("Save time = %f sec\n", save_finish - save_start);

	// Print out the output data if n is small.
	if (n <= printable_n)
	{
		printf("Output data:\n");
		for (int i = 0; i < n; i++)
			printf("(%f, %f)\n", h_pairs[i].first, h_pairs[i].second);
		printf("\n");
	}
#if 1
	// Show last 10 pairs.
	printf("Last 10 pairs:\n");
	for (int i = n - 11; i < n; i++)
	{
		printf("(%f, %f)\n", pairs[i].first, pairs[i].second);
	}
	printf("\n");
#endif
	return 0;
}

