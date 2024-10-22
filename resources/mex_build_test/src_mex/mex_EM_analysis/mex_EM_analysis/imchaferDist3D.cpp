//#include "mex.h"
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
using namespace std;
void calChaferDist3D(bool* mask1, bool* source,double* output_dist, size_t lenx, size_t leny, size_t lenz, size_t resx, size_t resy, size_t resz) {
	bool modif = true;
	double newVal;
	double valtemp;
	int nIter = 1;
	int iNew, jNew, kNew;
	int di1[13] = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, -1, -1 };
	int dj1[13] = { -1, -1, -1, 0, 0, 0, 1, 1,1,-1,-1,0,1 };
	int dk1[13] = { -1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,0, 0 };

	int di2[13] = { 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, 1,1 };
	int dj2[13] = { 1, 1,1, 0,0,0,-1,-1,-1,1,1,0,-1 };
	int dk2[13] = { 1,1,1,1,1,1,1,1,1,0,0,0,0 };
	double w1 = resx;
	double w2 = sqrt(2) * resx;
	double w3 = resz;
	double w4 = sqrt(resx * resx + resz * resz);
	double w5 = sqrt(resx * resx + resx * resx + resz * resz);
	double ws[13] = { w5, w4, w5, w4, w3, w4, w5,w4,w5,w2,w1,w1,w2 };
	while (modif)
	{
		modif = false;
		for (size_t k = 1; k < lenz - 1; k++)
		{
			for (size_t j = 1; j < leny - 1; j++)
			{
				for (size_t i = 1; i < lenx - 1; i++)
				{

					if (mask1[i + j * lenx + k * (lenx * leny)])
					{
						newVal = output_dist[i + j * lenx + k * (lenx * leny)];
						for (size_t v = 0; v < 13; v++)
						{
							iNew = i + di1[v];
							jNew = j + dj1[v];
							kNew = k + dk1[v];
							valtemp = output_dist[iNew + jNew * lenx + kNew * (lenx * leny)] + ws[v];
							newVal = min(newVal, valtemp);
						}

						if (newVal != output_dist[i + j * lenx + k * (lenx * leny)])
						{
							modif = true;
							output_dist[i + j * lenx + k * (lenx * leny)] = newVal;
						}
					}

				}
			}
		}
		if (modif == false && nIter != 1)
		{
			break;
		}
		modif = false;
		for (size_t k = lenz - 2; k > 0; k--)
		{
			for (size_t j = leny - 2; j > 0; j--)
			{
				for (size_t i = lenx - 2; i > 0; i--)
				{

					if (mask1[i + j * lenx + k * (lenx * leny)])
					{
						newVal = output_dist[i + j * lenx + k * (lenx * leny)];
						for (size_t v = 0; v < 13; v++)
						{
							iNew = i + di2[v];
							jNew = j + dj2[v];
							kNew = k + dk2[v];
							valtemp = output_dist[iNew + jNew * lenx + kNew * (lenx * leny)] + ws[v];
							newVal = min(newVal, valtemp);
						}

						if (newVal != output_dist[i + j * lenx + k * (lenx * leny)])
						{
							modif = true;
							output_dist[i + j * lenx + k * (lenx * leny)] = newVal;
						}
					}

				}
			}
		}
		if (modif == false && nIter != 1)
		{
			break;
		}
		nIter = nIter + 1;
	}
}

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
    // inputs: (logical array) foregroundMask, (logical array) source of distance searching, lenx, leny, lenz, resx, resy ,resz
	size_t lenx, leny, lenz, resx, resy, resz;
	bool* mask1;
	bool* source;
	mask1 = mxGetLogicals(prhs[0]);
	source = mxGetLogicals(prhs[1]);
	lenx = mxGetScalar(prhs[2]);
	leny = mxGetScalar(prhs[3]);
	lenz = mxGetScalar(prhs[4]);
	resx = mxGetScalar(prhs[5]);
	resy = mxGetScalar(prhs[6]);
	resz = mxGetScalar(prhs[7]);
	mxArray* output_dist_mx = mxCreateDoubleMatrix((mwSize)lenx*leny*lenz, (mwSize)1, mxREAL);
	double* output_dist = mxGetPr(output_dist_mx);
	for (size_t i = 0; i < lenx * leny * lenz; i++)
	{
		if (source[i])
		{
			output_dist[i] = 0;

		}
		else
		{
			output_dist[i] = 1e10;
		}
	}
	calChaferDist3D(mask1, source, output_dist, lenx, leny, lenz, resx, resy, resz);
	plhs[0] = output_dist_mx;

}