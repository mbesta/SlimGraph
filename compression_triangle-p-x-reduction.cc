// Copyright (c) 2019 ETH-Zurich.
//                    All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "benchmark.h"
#include "command_line.h"

#include <cassert>

#include <omp.h>

#define LOCKS_FACTOR 10

double SG_fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void SG_p_1_reduction( NodeID u, NodeID v, NodeID w, double tr_stays, int64_t locks_factor,
  std::vector<int64_t>** edges_to_be_removed, omp_nest_lock_t* locks ) {

  double rnd = SG_fRand(0, 1);
  if( tr_stays < rnd ) {
    // triangle will be compressed
    // a random one of the three edges will be removed
    rnd = SG_fRand(0, 1);
    if( (0 <= rnd) && (rnd < (1/3)) ) {
      // remove edge (u, v)
      omp_set_nest_lock(&locks[u / locks_factor]);
      omp_set_nest_lock(&locks[v / locks_factor]);

      edges_to_be_removed[u]->push_back( v );
      edges_to_be_removed[v]->push_back( u );

      omp_unset_nest_lock(&locks[u / locks_factor]);
      omp_unset_nest_lock(&locks[u / locks_factor]);
    } else {
      if( ((1/3) <= rnd) && (rnd < (2/3)) ) {
        // remove edge (v, w)
        omp_set_nest_lock(&locks[v / locks_factor]);
        omp_set_nest_lock(&locks[w / locks_factor]);

        edges_to_be_removed[v]->push_back( w );
        edges_to_be_removed[w]->push_back( v );

        omp_unset_nest_lock(&locks[v / locks_factor]);
        omp_unset_nest_lock(&locks[w / locks_factor]);
      } else {
        // remove edge (u, w)
        omp_set_nest_lock(&locks[u / locks_factor]);
        omp_set_nest_lock(&locks[w / locks_factor]);

        edges_to_be_removed[u]->push_back( w );
        edges_to_be_removed[w]->push_back( u );

        omp_unset_nest_lock(&locks[u / locks_factor]);
        omp_unset_nest_lock(&locks[w / locks_factor]);
      }
    }
  }
}

int main(int argc, char* argv[]) {
  // seed random number generator
  srand((unsigned)time(0));

  CLApp cli(argc, argv, "triangle compression kernel (p-x-reduction)");
  if (!cli.ParseArgs()) {
    return -1;
  }

  std::string g_output_file_name = cli.output_filename();
  double tr_stays = cli.sampling_propability();

  assert( (tr_stays >= 0) && (tr_stays <= 1) );

  // build a weighted graph
  Builder b(cli);
  Graph g = b.MakeGraph();

  if (g.directed()) {
    std::cout << "Input graph is directed but the p-x-reduction kernel from the triangle compression requires an undirected graph." << std::endl;
    return -2;
  }

  // general idea for the triangle compression algorithm:
  //
  // first we find all triangles
  //
  // for the triangles: run compression method,
  // and add edges to be removed to a global data structure
  //
  // create a copy of the existing graph, and remove the edges marked for removal
  //
  // write out the resulting graph
  //
  // performance-wise it would probably be better to mix the steps above,
  // but for clarity of the code, so it easier to understand, we do each step separately

  // find all triangles:
  // collect triangles in a global array so we can do each step separately
  // 
  // alternative idea:
  // it would be also possible to save the triangles locally in each thread,
  // and then process them in the same thread
  std::vector<NodeID>* triangles = new std::vector<NodeID>();

#pragma omp parallel for schedule(dynamic, 64)
  for (NodeID u = 0; u < g.num_nodes(); u++) {
    for (NodeID v : g.out_neigh(u)) {
      if (v > u) {
        // graph is undirected, so we want to avoid to look at each edge more than once
        break;
      }
      auto it = g.out_neigh(u).begin();
      for (NodeID w : g.out_neigh(v)) {
        if (w > v) {
          // graph is undirected, so we want to avoid to look at each edge more than once
          break;
        }
        // this works because we assume that the neighborhood vertices are sorted
        while (it != g.out_neigh(u).end() && (*it < w)) {
          it++;
        }
        if (w == *it) {
          // found a triangle
#pragma omp critical
          {
            triangles->push_back( u );
            triangles->push_back( v );
            triangles->push_back( w );
          }
        }
      }
    }
  }

  // init data structures for the edges that will be marked for removal
  std::vector<int64_t>** edges_to_be_removed = new std::vector<int64_t>*[g.num_nodes()];
#pragma omp parallel for
  for(int64_t v = 0; v < g.num_nodes(); ++v) {
    edges_to_be_removed[v] = new std::vector<int64_t>();
  }

  int locks_factor = LOCKS_FACTOR;
  int64_t locks_nr = (g.num_nodes() / locks_factor) + 1;
  omp_nest_lock_t* locks = new omp_nest_lock_t[locks_nr]();

#pragma omp parallel for
  for(int64_t i = 0; i < locks_nr; ++i) {
    omp_init_nest_lock(&locks[i]);
  }

  // go in parallel through all triangles and mark some of its edges for removal 
#pragma omp parallel for schedule(dynamic, 64)
  for(size_t i = 0; i < triangles->size(); i=i+3) {
    SG_p_1_reduction( (*triangles)[i], (*triangles)[i+1], (*triangles)[i+2], tr_stays, locks_factor, edges_to_be_removed, locks );
  }

#pragma omp parallel for
  for(int64_t i = 0; i < locks_nr; ++i) {
    omp_destroy_nest_lock(&locks[i]);
  }
  delete [] locks;

  delete triangles;

  // go in parallel through the uncompressed graph and only copy those edges that will stay
  std::vector<int64_t>** sampled_edges = new std::vector<int64_t>*[g.num_nodes()];
#pragma omp parallel for
  for(int64_t v = 0; v < g.num_nodes(); v++) {
    sampled_edges[v] = new std::vector<int64_t>();
  }

#pragma omp parallel for schedule(dynamic, 64)
  for (NodeID u = 0; u < g.num_nodes(); u++) {
    // first we have to sort the edges marked for removal
    //
    // it is also possible that an edge has been marked for removal more than once
    // (if it is in more than one triangle), so we make sure, that every edge is
    // only be found once in the datastructure
    std::sort(edges_to_be_removed[u]->begin(), edges_to_be_removed[u]->end());
    edges_to_be_removed[u]->erase( unique(edges_to_be_removed[u]->begin(), edges_to_be_removed[u]->end()), edges_to_be_removed[u]->end());

    // not sure that for( NodeID v : g.out_neigh(u)) gives me the elements in a sorted order
    auto it_neigh = g.out_neigh(u).begin();
    auto it_removal = edges_to_be_removed[u]->begin();
    while( it_neigh != g.out_neigh(u).end() ) {
      if( (it_removal != edges_to_be_removed[u]->end()) && (*it_neigh == *it_removal) ) {
        // found an edge that is marked for removal
        it_removal++;
      } else {
        sampled_edges[u]->push_back( *it_neigh );
      }
      it_neigh++;
    }
  }

#pragma omp parallel for
  for(int64_t i = 0; i < g.num_nodes(); ++i) {
    delete edges_to_be_removed[i];
  }
  delete [] edges_to_be_removed;

  // write out the compressed graph
  std::ofstream compressed_graph_file(g_output_file_name);

  for (NodeID i = 0; i < g.num_nodes(); i++) {
    for (size_t j=0 ; j < sampled_edges[i]->size() ; j++) {
      compressed_graph_file << i << " " << (*(sampled_edges[i]))[j] << std::endl;
    }
  }
 
#pragma omp parallel for
  for(int64_t i = 0; i < g.num_nodes(); ++i) {
    delete sampled_edges[i];
  }
  delete [] sampled_edges;

  return 0;
}
