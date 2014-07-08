#include "data.h"

#include <stdio.h>

// fields that hold the solution
double *x_new = NULL, *x_old = NULL; // 2d
double *bndN = NULL, *bndE = NULL, *bndS = NULL, *bndW = NULL;

struct discretization_t options;

