Slim Graph is a framework for practical lossy graph compression.

Slim Graph was written by Maciej Besta, Simon Weber, Lukas Gianinazzi,
Robert Gerstenberger, Andrey Ivanov, Yishai Oltchik and Torsten Hoefler
and is copyrighted by ETH-Zurich (c) 2019. See LICENSE for details on
the license.


Building
========

Dependencies
------------
Slim Graph is based in part on the GAP Benchmark Suite (GAPBS) and will
require an unpacked folder of their sources. Please specify the path to
the main GAPBS directory on your system in the file Makefile.inc
(variable GAPBS).

Two files in the Slim Graph folder are directly copied (benchmark.h and
builder.h) from GAPBS and command_line.h was also copied from their
code, but was also modified.

For more information on the GAP Benchmark Suite, see their homepage
http://gap.cs.berkeley.edu/benchmark.html or their paper (S. Beamer, K.
Asanović, and D. Patterson: The GAP Benchmark Suite (2015)). The source
code can be found on their github page: https://github.com/sbeamer/gapbs

Includes
--------
Set the variable GAPBS in Makefile.inc, so it points to the main
directory of your GAP Benchmark Suite installation.

Build parameter
---------------
Additionally you can set the C++ compiler in the Makefile.inc. While
developing Slim Graph, we used gcc, but other compiler should work as
well.


Using Slim Graph
================

Right now, Slim Graph offers different kernels for lossy graph
compression. Each executable will take an input file, apply the
compression kernel, and write the resulting graph to an output file in
the edgelist (.el) format. Slim Graph supports all file formats for the
input file, that the GAP Benchmark Suite (GAPBS) supports as well. Each
compression kernel is located in a different file and comes with a
different set of command line parameters. At the moment it is only
possible to apply one compression kernel at the time, but different
compression kernels can be applied one after another. Slim Graph doesn't
support weighted graphs at the moment.

Edgelist files (.el) have a very basic format: it contains two columns
with are divided by a space character. The first column contains the UID
of the origin vertex (a node in GAPBS language) and the second column
the UID of the target vertex. Each row of the file forms an edge.

The output file is specified with -o for all executables. The input
graph can either be specified with the help of an input file (parameter
-f) or can be generated on the fly with the help of GAPBS (see their
documentation for further help). All possible parameters can be listed
with <executable> -h, but not all parameters are significant for each
executable.

To our knowledge, GAPBS doesn't necessarily recognizes if an input graph
from a file is undirected. In such a case it can help to supply the
parameter -s.

In the following, we will briefly describe each available compression
kernel and their input parameter. For further information, please
consult the Slim Graph paper.

Vertex Kernel
-------------
The vertex kernel (compression_vertex) will remove all vertices with a
degree of zero or one. Afterwards the vertices will be enumerated anew,
so that their UIDs are consecutive. This kernel reduces the number of
vertices (as well as edges) in contrast to the other kernels, that only
reduce the amount of edges.
The vertex kernel takes no additional parameter.

compression_vertex -f <input file> -o <output file>

Edge Kernels
------------
Slim Graph provides two edge kernels: spectral sparsification and
uniform sampling.

Uniform sampling decides for each edge independently, whether that edge
stays in the compressed graph or not.
The user provides a probability treshold (-p) for that decision. -p
takes a floating point parameter between 0.0 and 1.0, where a higher
fraction means that more edges will remain in the compressed graph. The
default value for the probability treshold is 1.0, so no edge will be
removed.

compression_edge_uniform_sampling -f <input file> -o <output file> -p 0.5

Spectral sparsification also decides with a certain probability treshold
whether an edge will remain in the compressed graph or not. The treshold
is different for each edge, since it takes the degree of the incident
vertices into account. 
edge_stays = min( 1.0 , Υ / min( u.deg , v.deg ) )
where u is the origin vertex of the edge, and v the target vertex.

The fraction Y of remaining edges incident to each vertex can be
proportional to log(n) (Y = p * log(n), parameter -V log), where n is
the number of vertices in the input graph and p is a user specified
parameter (parameter -p with a floating point value between 0.0 and
1.0, with the default value 1.0). The other variant uses the average
degree as Y (Y = p * m/n, parameter -V avg), where m is the number of
edges in the input graph.

compression_edge_spectral_sparsify -f <input file> -o <output file> -p 0.5 -V [avg/log]

Triangle Kernels
----------------
Slim Graph provides one triangle kernel: p-1-reduction.

p-1-reduction decides for each triangle independently, whether that
triangle will be compressed for the compressed graph or not. If the
triangle will be compressed, then in case of p-1-reduction exactly one
edge will be removed. The edge to be removed of the possible three will
be decided at random.
p-1-reduction only works with undirected graphs.
The user provides a probability treshold (-p) for the general
compression decision of the triangle. -p takes a floating point
parameter between 0.0 and 1.0, where a higher fraction means that more
triangles will remain in the compressed graph. The default value for the
probability treshold is 1.0, so no triangle will be removed.

compression_triangle-p-x-reduction -f <input file> -o <output file> -p 0.5 [-s]


Questions
=========

If you have additional questions, feel free to contact us under the
following e-mail address:
slimgraph@spcl.inf.ethz.ch


Citation
========

Any published work which uses this software should include the
following citation:

----------------------------------------------------------------------
M. Besta, S. Weber, L. Gianinazzi, R. Gerstenberger, A. Ivanov,
Y. Oltchik and T. Hoefler: Slim Graph: Practical Lossy Graph
Compression for Approximate Graph Processing, Storage, and Analytics.
In: Proceedings of the International Conference on High Performance
Computing, Networking, Storage and Analysis (SC'19)
----------------------------------------------------------------------
