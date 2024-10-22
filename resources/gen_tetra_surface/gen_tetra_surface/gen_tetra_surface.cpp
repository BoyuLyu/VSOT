#include "C:\Program Files\MATLAB\R2021a\extern\include\mex.h"
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
int sub2ind(int lenx, int leny, int lenz, int idx, int idy, int idz) {
    int s;
    s = idx + (idy - 1) * lenx + (idz - 1) * lenx * leny;
    return s;
}

tuple <int, int, int> ind2sub(int lenx, int leny, int lenz, int id) {
    int idx, idy, idz;
    idz = id / (lenx * leny) + 1;
    idy = (id % (lenx * leny)) / lenx + 1;
    idx = (id % (lenx * leny)) % lenx;
    return make_tuple(idx, idy, idz);
}
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
// input: the id of the binary mask, xxshift(192*1), yyshift, zzshift, lenx, leny, lenz
// output: the generated surface of the object and all the surfaces of the tetrahedra
    double* coordinate = mxGetPr(prhs[0]);
    mwSize NElem = mxGetNumberOfElements(prhs[0]);
    double* xxshift = mxGetPr(prhs[1]);
    mwSize Nshift = mxGetNumberOfElements(prhs[1]);
    double* yyshift = mxGetPr(prhs[2]);
    double* zzshift = mxGetPr(prhs[3]);
    int lenx = mxGetScalar(prhs[4]);
    int leny = mxGetScalar(prhs[5]);
    int lenz = mxGetScalar(prhs[6]);
    int* idx = new int[NElem];
    int* idy = new int[NElem];
    int* idz = new int[NElem];
    int coordinatetmp;
    int* idxss = new int[Nshift* NElem]; // forming into a 2D array NElem x Nshift
    int* idyss = new int[Nshift* NElem];
    int* idzss = new int[Nshift* NElem];
    for (size_t i = 0; i < NElem; i++)
    {   
        coordinatetmp = int(coordinate[i]);
        tie(idx[i], idy[i], idz[i]) = ind2sub(lenx, leny, lenz, coordinatetmp);
        for (size_t j = 0; j < Nshift; j++)
        {
            idxss[i * Nshift + j] = idx[i] + xxshift[j];
            idyss[i * Nshift + j] = idy[i] + yyshift[j];
            idzss[i * Nshift + j] = idz[i] + zzshift[j];
        }
    }

























}