#define _USE_MATH_DEFINES
#include "stdio.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include "C:\Program Files\MATLAB\R2021a\extern\include\mex.h"
#include <algorithm>
#include <list>
#include <vector>
#include <queue>
#include <iterator>
#include <unordered_set>
#include <array>
using namespace std;

bool comp(vector<int>& rowA, vector<int>& rowB) {
    // support at most 4 columns
    int len = rowA.size();
    switch (len)
    {
        case 1:
            return(rowA[0] < rowB[0]);
            break;
        case 2:
            return(rowA[0] < rowB[0] || (rowA[0] == rowB[0] && rowA[1] < rowB[1]));
            break;

        case 3:
            return(rowA[0] < rowB[0] || (rowA[0] == rowB[0] && rowA[1] < rowB[1]) || (rowA[0] == rowB[0] && rowA[1] == rowB[1] && rowA[2] < rowB[2]));
            break;

        case 4:
            return(rowA[0] < rowB[0] || (rowA[0] == rowB[0] && rowA[1] < rowB[1]) || (rowA[0] == rowB[0] && rowA[1] == rowB[1] && rowA[2] < rowB[2]) || (rowA[0] == rowB[0] && rowA[1] == rowB[1] && rowA[2] == rowB[2] && rowA[3] < rowB[3]));
            break;
    }

}


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {

    double* inputMat;
    inputMat = mxGetPr(prhs[0]);
    int m = mxGetM(prhs[0]);
    int n = mxGetN(prhs[0]);
    //mexPrintf("\n%d rows %d col\n", m, n);
    mxArray* dtMex = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    double* dt1 = mxGetPr(dtMex);
    int newM;
    array <double, 4> rowx = {0,0,0,0};
    switch (n)
    {

    case 2:
        //sort each row
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n ; j++) {
                rowx[j] = inputMat[i + j * m];
            }
            sort(rowx.begin(), rowx.begin() + 2);
            for (size_t j = 0; j < n; j++) {
                dt1[i + j*m] = rowx[j];
            }

        }
        break;
    case 3:

        //sort each row
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n ; j++) {
                rowx[j] = inputMat[i + j * m];
            }
            sort(rowx.begin(), rowx.begin() + 3);
            for (size_t j = 0; j < n; j++) {
                dt1[i + j * m] = rowx[j];
            }

        }
        break;

    case 4: 
        //sort each row
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++) {
                rowx[j] = inputMat[i + j * m];
            }
            sort(rowx.begin(), rowx.begin() + 4);
            for (size_t j = 0; j < n; j++) {
                dt1[i + j * m] = rowx[j];
            }

        }
        break;
    }

    plhs[0] = dtMex;
    // //sort rows based on the column
    //switch (m)
    //{
    //    case 3:
    //        // the last column is the index of the triangle, first two column represent the edge






    //        break;
    //    case 4:
    //        // the last column is the index of the tetrahedron, the first three represents one of the triangles in it

    //        break;
    //}

    //vector<vector<int>> mat(m, vector<int>(n, 0));
    //for (size_t i = 0; i < m; i++)
    //{
    //    vector <int> tmp(n,0);
    //    for (size_t j = 0; j < n; j++)
    //    {
    //        tmp[j] = inputMat[i + j * m];
    //    }
    //    //mexPrintf("\n%d length\n", tmp.size(), tmp[1]);
    //    sort(tmp.begin(), tmp.end());
    //    mat[i] = tmp;
    //}

    //sort(mat.begin(), mat.end(), &comp);
    //vector< vector<int> >::iterator new_end = unique(mat.begin(), mat.end());
    //newM = distance(mat.begin(), new_end);
    //mat.resize(distance(mat.begin(), new_end));
    //mxArray* dt2Mex = mxCreateDoubleMatrix((mwSize)newM, (mwSize)n, mxREAL);
    //double* dt2 = mxGetPr(dt2Mex);
    //int count = 0;
    //for (auto it = mat.begin(); it != mat.end(); ++it)
    //{
    //    vector<int> tmp = *it;
    //    for (size_t i = 0; i < n; i++)
    //    {
    //        dt2[count + i * newM] = tmp[i];
    //    }
    //    count++;


    //}
    //plhs[0] = dt2Mex;

}



