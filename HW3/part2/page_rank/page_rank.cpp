#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>
#include <vector>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{

  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  for (int i = 0; i < numNodes; ++i)
  {
    solution[i] = equal_prob; //initialize
  }

  /*
     For PP students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */

  // allocate
  double *score_new = new double[numNodes];
  double *score_out = new double[numNodes];
  int *num_out = new int[numNodes];
  int *num_in = new int[numNodes];
  std::vector<int> end_vertex;

#pragma omp parallel for schedule(guided)
  for (int i = 0; i < numNodes; ++i)
  {
    num_out[i] = outgoing_size(g, i);
    if (num_out[i] == 0){ // if the vertax has no out edge
#pragma omp critical
      end_vertex.push_back(i);
    }
    num_in[i] = incoming_size(g, i);
  }

  double global_diff = convergence;
  double add_up = (1.0 - damping) / numNodes;

  while (global_diff >= convergence){ // loop until converge
    // compute out value
#pragma omp parallel for
    for (int i = 0; i < numNodes; ++i){
      if (num_out[i] != 0) score_out[i] = solution[i] / num_out[i];
    }

    // updating new score
#pragma omp parallel for schedule(guided)
    for (int i = 0; i < numNodes; ++i){
      score_new[i] = 0;
      for (int j = 0; j < num_in[i]; ++j){
        score_new[i] += score_out[g->incoming_edges[g->incoming_starts[i] + j]];
      }
      score_new[i] = score_new[i] * damping + add_up;
    }

    // add score from no out edge vertax
    double s = 0;
#pragma omp parallel for reduction(+ \
                                   : s)
    for (int i = 0; i < end_vertex.size(); ++i){
      s += damping * solution[end_vertex[i]] / numNodes;
    }

#pragma omp parallel for
    for (int i = 0; i < numNodes; ++i){
      score_new[i] += s;
    }

    // compute global difference
    global_diff = 0;
#pragma omp parallel for reduction(+ \
                                   : global_diff)
    for (int i = 0; i < numNodes; ++i){
      global_diff += abs(solution[i] - score_new[i]);
    }

    for (int i = 0; i < numNodes; ++i){
      solution[i] = score_new[i];
    }
  }

  // free 
  delete[] score_new;
  delete[] score_out;
  delete[] num_out;
  delete[] num_in;

}
