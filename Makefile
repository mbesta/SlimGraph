include Makefile.inc

CXX_FLAGS += -std=c++11 -O3 -Wall -Wextra -fopenmp -I$(GAPBS)/src

NEXT_RELEASE=0.1

EXE = \
	compression_edge_spectral_sparsify \
	compression_edge_uniform_sampling \
	compression_triangle-p-x-reduction \
	compression_vertex

.PHONY: all
all: $(EXE)

.PHONY: clean
clean:
	rm -f $(EXE)

compression_edge_spectral_sparsify: compression_edge_spectral_sparsify.cc command_line.h
	$(CXX) $(CXX_FLAGS) $< -o $@

compression_edge_uniform_sampling: compression_edge_uniform_sampling.cc command_line.h
	$(CXX) $(CXX_FLAGS) $< -o $@

compression_triangle-p-x-reduction: compression_triangle-p-x-reduction.cc 
	$(CXX) $(CXX_FLAGS) $< -o $@

compression_vertex: compression_vertex.cc command_line.h
	$(CXX) $(CXX_FLAGS) $< -o $@
