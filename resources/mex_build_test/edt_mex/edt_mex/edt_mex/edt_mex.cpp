#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <algorithm>
#include <list>
#include <vector>
#include <queue>
#include <iterator>
#include <unordered_set>
#include "edt.hpp"
using namespace std;
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
    bool* labels3d = mxGetLogicals(prhs[0]);
    int lenx, leny, lenz, resx, resy, resz;
	lenx = mxGetScalar(prhs[1]);
	leny = mxGetScalar(prhs[2]);
	lenz = mxGetScalar(prhs[3]);
	resx = mxGetScalar(prhs[4]);
	resy = mxGetScalar(prhs[5]);
	resz = mxGetScalar(prhs[6]);
	mxArray* dtMex = mxCreateDoubleMatrix((mwSize)lenx * leny * lenz, (mwSize)1, mxREAL);
	double* dt1 = mxGetPr(dtMex);
	float* dt = edt::edt<bool>(labels3d,
		/*sx=*/lenx, /*sy=*/leny, /*sz=*/lenz,
		/*wx=*/resx, /*wy=*/resy, /*wz=*/resz,
		/*black_border=*/true);
	for (size_t i = 0; i < lenx * leny * lenz; i++)
	{
		dt1[i] = double(dt[i]);
	}
	plhs[0] = dtMex;

}
