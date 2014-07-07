#include <cstdio>
#include <cstdlib>
#include <time.h>

#ifndef CPU
#pragma offload_attribute(push,target(mic))
#endif

#include <algorithm>
#include <tbb/parallel_invoke.h>
#include <utility>

// merge sequences [xs,xe) and [ys,ye) to output [zs,(xe-xs)+(ye-ys)
template<typename T>
static void parallel_merge(const T* xs, const T* xe, const T* ys, const T* ye, T* zs)
{
	const size_t MERGE_CUT_OFF = 2000;
	if (xe - xs + ye - ys <= MERGE_CUT_OFF)
	{
		std::merge(xs, xe, ys, ye, zs);
	}
	else
	{
		const T *xm, *ym;
		if (xe - xs < ye - ys)
		{
			ym = ys + (ye - ys) / 2;
			xm = std::upper_bound(xs, xe, *ym);
		}
		else
		{
			xm = xs + (xe - xs) / 2;
			ym = std::lower_bound(ys, ye, *xm);
		}
		T* zm = zs + (xm - xs) + (ym - ys);
		tbb::parallel_invoke(
			[=] { parallel_merge(xs, xm, ys, ym, zs); },
			[=] { parallel_merge(xm, xe, ym, ye, zm); } );
    }
}

// sorts [xs,xe).  zs[0:xe-xs) is temporary buffer supplied by caller.
// result is in [xs,xe) if inplace==true, otherwise in zs[0:xe-xs)
template<typename T>
void parallel_merge_sort(T* xs, T* xe, T* zs, bool inplace)
{
	const size_t SORT_CUT_OFF = 500;
	if (xe - xs <= SORT_CUT_OFF)
	{
		std::stable_sort(xs, xe);
		if (!inplace)
			std::move(xs, xe, zs);
	}
	else
	{
		T* xm = xs + (xe - xs) / 2;
		T* zm = zs + (xm - xs);
		T* ze = zs + (xe - xs);
		
		tbb::parallel_invoke(
			[=] { parallel_merge_sort(xs, xm, zs, !inplace); },
			[=] { parallel_merge_sort(xm, xe, zm, !inplace); });
		
		if (inplace)
			parallel_merge(zs, zm, zm, ze, xs);
		else
			parallel_merge(xs, xm, xm, xe, zs);
	}
}

#ifndef CPU
#pragma offload_attribute(pop)
#endif

// Get the timer float.
static void get_time(double* ret)
{
	volatile struct timespec val;
	clock_gettime(CLOCK_REALTIME, (struct timespec*)&val);
	*ret = (double)0.000000001 * val.tv_nsec + val.tv_sec;
}

void usage(const char* filename)
{
	printf("Sort the random key-float data set of the given length by key.\n");
	printf("Usage: %s <n>\n", filename);
}

using namespace tbb;

// Memory alignment, for vectorization on MIC.
// 4096 should be best for memory transfers over PCI-E.
#define MEMALIGN 4096

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

	std::pair<float, float>* pairs; posix_memalign(
		(void**)&pairs, MEMALIGN, n * sizeof(std::pair<float, float>));

	// Generate the random keys and floats on the host.
	for (int i = 0; i < n; i++)
	{
		pairs[i].first = drand48();
		pairs[i].second = drand48();
	}

	// Print out the input data if n is small.
	if (n <= printable_n)
	{
		printf("Input data:\n");
		for (int i = 0; i < n; i++)
			printf("(%f, %f)\n", pairs[i].first, pairs[i].second);
		printf("\n");
	}

	// Strip initialization delay.
	double init_start, init_finish;
	get_time(&init_start);
#ifndef CPU
	#pragma offload target(mic) \
		nocopy(pairs:length(n) alloc_if(0) free_if(0))
	{ }
#endif
	get_time(&init_finish);

	// Transfer data to the device.
	double load_start, load_finish;
	get_time(&load_start);
#ifndef CPU
	#pragma offload target(mic) \
		in(pairs:length(n) alloc_if(1) free_if(0))
	{ }
#endif
	get_time(&load_finish);

	std::pair<float, float>* temp;
#ifndef CPU
	#pragma offload target(mic)
#endif
	{	
		temp = new std::pair<float, float>[n];
	}

	double sort_start, sort_finish;
	get_time(&sort_start);
#ifndef CPU
	#pragma offload target(mic) \
		nocopy(pairs:length(n) alloc_if(0) free_if(0))
#endif
	{
		parallel_merge_sort(pairs, pairs + n, temp, true);
	}
	get_time(&sort_finish);

#ifndef CPU
	#pragma offload target(mic)
#endif
	{
		delete[] temp;
	}

	// Transfer data back to host.
	double save_start, save_finish;
	get_time(&save_start);
#ifndef CPU
	#pragma offload target(mic) \
		out(pairs:length(n) alloc_if(0) free_if(1))
	{ }
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
			printf("(%f, %f)\n", pairs[i].first, pairs[i].second);
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

