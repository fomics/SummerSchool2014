#include <dlfcn.h>
#include <stdio.h>

extern "C" int cudaLaunch_(const char* name)
{
	void* handle = dlopen(NULL, RTLD_NOW);
	if (!handle)
	{
		fprintf(stderr, "%s\n", dlerror());
		exit(-1);
	}
	void* func = dlsym(handle, name);
	if (!func)
	{
		fprintf(stderr, "%s\n", dlerror());
		exit(-2);
	}

	return cudaLaunch(func);
}
