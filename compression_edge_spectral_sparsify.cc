// Copyright (c) 2019 ETH-Zurich.
//                    All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "benchmark.h"
#include "command_line.h"

#include <algorithm>
#include <cassert>

#include <omp.h>

#define LOCKS_FACTOR 10

double SG_fRand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

void SG_spectral_sparsify_directed( NodeID origin, NodeID target, double Y, const Graph &g, std::vector<int64_t>** sampled_edges ) {
  double rnd = SG_fRand(0, 1);
  double edge_stays = std::min(1.0, Y / std::min(g.out_degree(origin), g.in_degree(target)));
  if( edge_stays >= rnd ) {
    sampled_edges[origin]->push_back(target);
  }
}

void SG_spectral_sparsify_undirected( NodeID u, NodeID v, double Y, int64_t locks_factor,
  const Graph &g, std::vector<int64_t>** sampled_edges, omp_nest_lock_t* locks ) {
  
  double rnd = SG_fRand(0, 1);
  double edge_stays = std::min(1.0, Y / std::min(g.out_degree(u), g.in_degree(v)));
  if( edge_stays >= rnd ) {
    omp_set_nest_lock(&locks[u / locks_factor]);
    omp_set_nest_lock(&locks[v / locks_factor]);

    sampled_edges[u]->push_back(v);
    sampled_edges[v]->push_back(u);

    omp_unset_nest_lock(&locks[u / locks_factor]);
    omp_unset_nest_lock(&locks[v / locks_factor]);
  }
}

int main(int argc, char* argv[]) {
  // seed random number generator
  srand((unsigned)time(0));

  CLApp cli(argc, argv, "single-edge compression kernel (random uniform)");
  if (!cli.ParseArgs()) {
    return -1;
  }

  std::string g_output_file_name = cli.output_filename();
  double p = cli.sampling_propability();
  std::string spectral_variant = cli.spectral_variant();

  assert( (p >= 0) && (p <= 1) );

  // build graph
  Builder b(cli);
  Graph g = b.MakeGraph();

  int64_t n = g.num_nodes();

  // determine connectivity spectral parameter Y
  double Y;
  if( spectral_variant == "log" ) {
    Y = p * log2((double)n);
  } else {
    if( spectral_variant == "avg" ) {
      double avg_deg = 0;
#pragma omp parallel for reduction(+ : avg_deg)
      for(int64_t v = 0; v < n; ++v) {
        avg_deg += g.out_degree(v);
      }
      avg_deg /= (double)n;
      Y = p * avg_deg;
    } else {
      std::cerr << "Parameter v needs to have one of two values: 'avg' or 'log'." << std::endl;
      return -1;
    }
  }

  std::vector<int64_t>** sampled_edges = new std::vector<int64_t>*[n];
#pragma omp parallel for
  for(int64_t v = 0; v < n; ++v) {
    sampled_edges[v] = new std::vector<int64_t>();
  }

  if(g.directed()) {
#pragma omp parallel for schedule(dynamic, 64)
    for (NodeID u = 0; u < n; u++) {
      for (NodeID v : g.out_neigh(u)) {
        SG_spectral_sparsify_directed( u, v, Y, g, sampled_edges );
      }
    }
  } else {
    // graph is undirected
    int locks_factor = LOCKS_FACTOR;
    int64_t locks_nr = (n / locks_factor) + 1;
    omp_nest_lock_t* locks = new omp_nest_lock_t[locks_nr]();

#pragma omp parallel for
    for(int64_t i = 0; i < locks_nr; ++i) {
      omp_init_nest_lock(&locks[i]);
    }

#pragma omp parallel for schedule(dynamic, 64)
    for (NodeID u = 0; u < n; u++) {
      for (NodeID v : g.out_neigh(u)) {
        if (v < u) {
          SG_spectral_sparsify_undirected( u, v, Y, locks_factor, g, sampled_edges, locks );
        }
      }
    }

#pragma omp parallel for
    for(int64_t i = 0; i < locks_nr; ++i) {
      omp_destroy_nest_lock(&locks[i]);
    }
    delete [] locks;

  }

// write out the edges sorted: first by origin, then by target
#pragma omp parallel for
  for (int64_t i = 0; i < n; i++) {
    std::sort(sampled_edges[i]->begin(), sampled_edges[i]->end());
    // TODO: Why do we need that?
    // sampled_edges[i]->erase( unique(sampled_edges[i]->begin(), sampled_edges[i]->end() ), sampled_edges[i]->end() );
  }

// write out the compressed graph
  std::ofstream compressed_graph_file(g_output_file_name);
  
  for (NodeID i = 0; i < n; i++) {
    for (size_t j=0 ; j < sampled_edges[i]->size() ; j++) {
      compressed_graph_file << i << " " << (*(sampled_edges[i]))[j] << std::endl;
    }
  }

  compressed_graph_file.close();

#pragma omp parallel for
  for(int64_t i = 0; i < g.num_nodes(); ++i) {
    delete sampled_edges[i];
  }
  delete [] sampled_edges;

  return 0;
}
