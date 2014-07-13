void sincos_(int* nx, int* ny, int* nz);

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("Usage: %s <nx> <ny> <nz>\n", argv[0]);
		return 0;
	}

	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);

	sincos_(&nx, &ny, &nz);

	return 0;
}
