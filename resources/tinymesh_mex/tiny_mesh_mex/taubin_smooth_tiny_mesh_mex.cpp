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
#include "tinymesh.h"

namespace tms = tinymesh;



void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {

    double* vertices_mat = mxGetPr(prhs[0]);
    int num_vex = mxGetM(prhs[0]);
    double* indices_mat = mxGetPr(prhs[1]);
    int num_face = mxGetM(prhs[1]);
    int num_out_vex;
    int num_out_face;
    std::vector<Vec3> vertices;
    std::vector<uint32_T> indices;
    double shrink_smooth = mxGetScalar(prhs[2]);
    double inflate_smooth = mxGetScalar(prhs[3]);  
    double iters = mxGetScalar(prhs[4]);
    for (size_t i = 0; i < num_vex; i++)
    {
        Vec3 vertices_row(vertices_mat[i], vertices_mat[i + num_vex], vertices_mat[i + num_vex * 2], 0.0);
        vertices.push_back(vertices_row);

    }
    for (size_t i = 0; i < num_face; i++)
    {
        indices.push_back((uint32_T)indices_mat[i] -1 );
        indices.push_back((uint32_T)indices_mat[i + num_face]- 1);
        indices.push_back((uint32_T)indices_mat[i + 2*num_face] - 1);

    }

    const std::vector<Vec3> vertices_const = vertices;
    const std::vector<uint32_T> indices_const = indices;
    tms::Mesh mesh(vertices_const, indices_const);
    tms::smoothTaubin(mesh, shrink_smooth, inflate_smooth, (int)iters);
    const std::vector<Vec3> outvertices = mesh.getVertices();
    const std::vector<uint32_T> outindices = mesh.getVertexIndices();
    num_out_vex = outvertices.size();
    num_out_face = outindices.size();
    mxArray* out_vertices = mxCreateDoubleMatrix((mwSize) num_out_vex, 3, mxREAL);
    double* out_vertices_mat = mxGetPr(out_vertices);
    mxArray* out_faces = mxCreateDoubleMatrix((mwSize)num_out_face, 1, mxREAL);
    double* out_face_mat = mxGetPr(out_faces);
    for (size_t i = 0; i < num_out_vex; i++)
    {
        Vec3 tmp = outvertices[i];
        out_vertices_mat[i] = tmp[0];
        out_vertices_mat[i + num_out_vex] = tmp[1];
        out_vertices_mat[i + num_out_vex*2] = tmp[2];
    }
    for (size_t j = 0; j < num_out_face; j++)
    {
        out_face_mat[j] = outindices[j];
    }
    plhs[0] = out_vertices;
    plhs[1] = out_faces;













}