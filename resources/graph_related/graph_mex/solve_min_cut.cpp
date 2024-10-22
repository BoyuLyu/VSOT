#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "ibfs/ibfs.h"
#include "Graph.h"
#include <list>
#include "ibfs/instances.inc"

/*ibfs find min-cut*/
typedef double EnergyType;
mxClassID MATLAB_ENERGYTERM_TYPE = mxDOUBLE_CLASS;

typedef double EnergyTermType;
mxClassID MATLAB_ENERGY_TYPE = mxDOUBLE_CLASS;

typedef double LabelType;
mxClassID MATLAB_LABEL_TYPE = mxDOUBLE_CLASS;

typedef IBFSGraph<EnergyTermType, EnergyTermType, EnergyType> GraphType;

double round_mex(double a);
int isInteger(double a);

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif

/*min cut based on ibfs modified from
 *
 *Anton Osokin, (firstname.lastname@gmail.com)
 * https://github.com/aosokin/graphCutMex_IBFS
 ***********************************/

double round_mex(double a)
{
    return floor(a + 0.5);
}

int isInteger(double a)
{
    return (abs(a - round_mex(a)) < 1e-6);
}

void ibfsMinCut(mxArray* uInPtr, mxArray* pInPtr, mxArray** cOutPtr, mxArray** lOutPtr) {
    /*
    * uInPtr: unary
    * terminalWeights=[
    *16,0;
    *13,0;
    *0,20;
    *0,4
    *   ];
    * pInPtr: pairwise
    * edgeWeights=[
    * 1,2,10.9,4.5;
    * 1,3,12.1,-1.5;
    * 2,3,-1.5,9.3;
    * 2,3,14.2,0;
    * 1,4,0,7.5
    * ];
    * cOutPtr: cut
    * lOutPtr: label
    */
    // 
    // 
    // 
    //node number
    int numNodes;

    // get unary potentials
    MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "graphCutMex: The second paramater is not 2-dimensional");
    MATLAB_ASSERT(mxGetClassID(uInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Unary potentials are of wrong type");
    MATLAB_ASSERT(mxIsComplex(uInPtr) == false, "graphCutMex: Unary potentials should not be complex");

    numNodes = mxGetM(uInPtr);

    MATLAB_ASSERT(numNodes >= 1, "graphCutMex: The number of nodes is not positive");
    MATLAB_ASSERT(mxGetN(uInPtr) == 2, "graphCutMex: The second paramater is not of size #nodes x 2");

    EnergyTermType* termW = (EnergyTermType*)mxGetData(uInPtr);

    //get pairwise potentials
    MATLAB_ASSERT(mxGetNumberOfDimensions(pInPtr) == 2, "graphCutMex: The third paramater is not 2-dimensional");

    mwSize numEdges = mxGetM(pInPtr);

    MATLAB_ASSERT(mxGetN(pInPtr) == 4, "graphCutMex: The third paramater is not of size #edges x 4");
    MATLAB_ASSERT(mxGetClassID(pInPtr) == MATLAB_ENERGYTERM_TYPE, "graphCutMex: Pairwise potentials are of wrong type");

    EnergyTermType* edges = (EnergyTermType*)mxGetData(pInPtr);
    for (int i = 0; i < numEdges; ++i)
    {
        MATLAB_ASSERT(1 <= round_mex(edges[i]) && round_mex(edges[i]) <= numNodes, "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(isInteger(edges[i]), "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(1 <= round_mex(edges[i + numEdges]) && round_mex(edges[i + numEdges]) <= numNodes, "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(isInteger(edges[i + numEdges]), "graphCutMex: error in pairwise terms array");
        MATLAB_ASSERT(edges[i + 2 * numEdges] + edges[i + 3 * numEdges] >= 0, "graphCutMex: error in pairwise terms array: nonsubmodular edge");
    }


    // start computing
    //prepare graph
    GraphType* g = new GraphType(numNodes, numEdges);
    for (int i = 0; i < numNodes; ++i)
    {
        g->add_node(1);
        g->add_tweights(i, termW[i], termW[numNodes + i]);
    }

    for (int i = 0; i < numEdges; ++i)
        if (edges[i] < 1 || edges[i] > numNodes || edges[numEdges + i] < 1 || edges[numEdges + i] > numNodes || edges[i] == edges[numEdges + i] || !isInteger(edges[i]) || !isInteger(edges[numEdges + i])) {
            mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge has invalid vertex numbers and therefore it is ignored");
        }
        else
            if (edges[2 * numEdges + i] + edges[3 * numEdges + i] < 0) {
                mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Some edge is non-submodular and therefore it is ignored");
            }
            else
            {
                if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] >= 0)
                    g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), edges[2 * numEdges + i], edges[3 * numEdges + i]);
                else
                    if (edges[2 * numEdges + i] <= 0 && edges[3 * numEdges + i] >= 0)
                    {
                        g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i] + edges[2 * numEdges + i]);
                        g->add_tweights((GraphType::node_id)round_mex(edges[i] - 1), 0, edges[2 * numEdges + i]);
                        g->add_tweights((GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, -edges[2 * numEdges + i]);
                    }
                    else
                        if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] <= 0)
                        {
                            g->add_edge((GraphType::node_id)round_mex(edges[i] - 1), (GraphType::node_id)round_mex(edges[numEdges + i] - 1), edges[3 * numEdges + i] + edges[2 * numEdges + i], 0);
                            g->add_tweights((GraphType::node_id)round_mex(edges[i] - 1), 0, -edges[3 * numEdges + i]);
                            g->add_tweights((GraphType::node_id)round_mex(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i]);
                        }
                        else
                            mexWarnMsgIdAndTxt("graphCutMex:pairwisePotentials", "Something strange with an edge and therefore it is ignored");
            }

    //compute flow
    EnergyType flow = g->maxflow();

    //output minimum value
    if (cOutPtr != NULL) {
        *cOutPtr = mxCreateNumericMatrix(1, 1, MATLAB_ENERGY_TYPE, mxREAL);
        *(EnergyType*)mxGetData(*cOutPtr) = (EnergyType)flow;
    }

    //output minimum cut
    if (lOutPtr != NULL) {
        *lOutPtr = mxCreateNumericMatrix(numNodes, 1, MATLAB_LABEL_TYPE, mxREAL);
        LabelType* segment = (LabelType*)mxGetData(*lOutPtr);
        for (int i = 0; i < numNodes; i++)
            segment[i] = g->what_segment(i);
    }

    delete g;




}
/*Use graph to find the connected component using Graph.h*/
// Graph class represents a undirected graph
// using adjacency list representation


// Method to print connected components in an
// undirected graph
void Graph::connectedComponents()
{
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    int curLabel = 0;
    for (int v = 0; v < V; v++)
        visited[v] = false;

    for (int v = 0; v < V; v++) {
        if (visited[v] == false) {
            // print all reachable vertices
            // from v
            DFSUtil(v, visited, labelsOutput, curLabel);
            ++curLabel;
        }
    }
    delete[] visited;

}

void Graph::DFSUtil(int v, bool visited[], int labels[], int curLabel)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    labels[v] = curLabel;
    // Recur for all the vertices
    // adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i) {
        if (!visited[*i]) {
            DFSUtil(*i, visited, labels, curLabel);
        }
    }
}

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    labelsOutput = new int[V];
}

Graph::~Graph() { delete[] adj; }

// method to add an undirected edge
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {
    int numNode = mxGetM(prhs[0]);

    mxArray* labelMatmx = mxCreateDoubleMatrix(numNode, 1, mxREAL);
    mxArray* cutMatmx = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray* uInPtr = (mxArray*)prhs[0];
    mxArray* pInPtr = (mxArray*)prhs[1];
    ibfsMinCut(uInPtr, pInPtr, &cutMatmx, &labelMatmx);
    plhs[0] = labelMatmx;
}