// Copyright (c) 2019 ETH-Zurich.
//                    All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "benchmark.h"
#include "command_line.h"

#include <omp.h>

void SG_low_degree( NodeID u, const Graph &g, bool keep_nodes[] ) {
  if( (g.in_degree( u ) + g.out_degree( u )) < 2 ) { // equivalent to (v.deg == 0) || (v.deg == 1)
    keep_nodes[u] = false;
  } else {
    keep_nodes[u] = true;
  }
}

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "single-vertex compression kernel");
  if (!cli.ParseArgs()) {
    return -1;
  }

  std::string g_output_file_name = cli.output_filename();

  // build graph
  Builder b(cli);
  Graph g = b.MakeGraph();

  int64_t n = g.num_nodes();

  bool* keep_nodes = new bool[n];

#pragma omp parallel for
  for (NodeID i = 0; i < n; i++) {
    SG_low_degree( i, g, keep_nodes );
  }

  // we need to reassign the IDs of the vertices
  size_t num_removed = 0;
  NodeID* newNodeID = new NodeID[n];
  for (NodeID i = 0; i < n; i++) {
    if( keep_nodes[i] ) {
      newNodeID[i] = i - num_removed;
    } else {
      num_removed++;
    }
  }

  // write out the compressed graph
  std::ofstream compressed_graph_file(g_output_file_name);

  for (NodeID i = 0; i < n; i++) {
    if( keep_nodes[i] ) {
      for (NodeID j : g.out_neigh(i)) {
        if( keep_nodes[j] ) {
          compressed_graph_file << newNodeID[i] << " " << newNodeID[j] << std::endl;
        }
      }
    }
  }

  compressed_graph_file.close();

  delete [] keep_nodes;
  delete [] newNodeID;

  return 0;
}
