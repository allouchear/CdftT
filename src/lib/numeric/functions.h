#ifndef _CDFTT_FUNCTIONS_H_INCLUDED
#define _CDFTT_FUNCTIONS_H_INCLUDED

using namespace std;
#include <numeric/Grid.h>
#include <cmath>
#include <vector>

double coulomb(vector<double>, double, double, double, const Grid& );

double laplacian(vector<double>, double, double, double, const Grid&);

double grad(vector<double>, double, double, double, const Grid&);
#endif
