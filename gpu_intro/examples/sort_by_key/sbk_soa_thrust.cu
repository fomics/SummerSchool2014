#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

// Get the timer value.
void get_time(double* ret)
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
	host_vector<float> h_keys(n);
	host_vector<float> h_vals(n);
	for (int i = 0; i < n; i++)
	{
		h_keys[i] = drand48();
		h_vals[i] = drand48();
	}

	// Print out the input data if n is small.
	if (n <= printable_n)
	{
		printf("Input data:\n");
		for (int i = 0; i < n; i++)
			printf("(%f, %f)\n", h_keys[i], h_vals[i]);
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
	device_vector<float> d_keys = h_keys;
	device_vector<float> d_vals = h_vals;
#endif
	get_time(&load_finish);

	double sort_start, sort_finish;
	get_time(&sort_start);
#ifdef CPU
	sort_by_key(h_keys.begin(), h_keys.end(), h_vals.begin());
#else
	sort_by_key(d_keys.begin(), d_keys.end(), d_vals.begin());
	cudaDeviceSynchronize();
#endif
	get_time(&sort_finish);

	// Transfer data back to host.
	double save_start, save_finish;
	get_time(&save_start);
#ifndef CPU
	copy(d_keys.begin(), d_keys.end(), h_keys.begin());
	copy(d_vals.begin(), d_vals.end(), h_vals.begin());
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
			printf("(%f, %f)\n", h_keys[i], h_vals[i]);
		printf("\n");
	}
#if 1
	// Show last 10 pairs.
	printf("Last 10 pairs:\n");
	for (int i = n - 11; i < n; i++)
	{
		printf("(%f, %f)\n", h_keys[i], h_vals[i]);
	}
	printf("\n");
#endif
	return 0;
}

