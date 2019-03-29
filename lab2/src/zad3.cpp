#include "matrix.h"
#include "circuit.h"


int main(int argc, char **argv) {
    
    if (argc != 2)
		fprintf(stderr, "Usage: %s <connections file>\n", argv[0]);
    
    
    struct graph graph_inst;    
    load_graph(argv[1], &graph_inst);
	
	
	const int N = graph_inst.node_count;
	
	double A[N*N];
	int U[N];
	double I[N];
	
	for(int y=0; y<N; y++) {
		for(int x=0; x<N ; x++)
			A[x + N*y] = 0;
		I[y] = 0;
		U[y] = y;
	}
		
	prepare_equations(A, I, N, &graph_inst);	
	solve_matrix(A, U, I, N);
	show_matricies(A, U, I, 1, N);
	
	update_currents(U, I, &graph_inst);
	update_graph(&graph_inst, argv[1]);
	
	return 0;
}
