#include "circuit.h"


void update_currents(int U[], double I[], struct graph *graph_ptr) {
	// store results of solving matricies in structure
	for(int i=0; i<graph_ptr->node_count; i++)
		graph_ptr->nodes[U[i]]->U = I[i];
	
	//update currents in edges and store the biggest one by module
	struct edge **edges = graph_ptr->edges;
	for(int i=0; i<graph_ptr->edge_count; i++) {
		edges[i]->I = (edges[i]->start->U -edges[i]->end->U) / edges[i]->R;
	}
	
	// add all currents going out of the source node
	struct node *source = graph_ptr->source;
	double max_I = 0.0;
	for(int i=0; i<source->edges.size(); i++) {
		max_I += fabs(source->edges[i]->I);
	}
	
	// dont want to end up with something very small in max_I, becouse it's used for scaling
	graph_ptr->max_I = max_I > 0.005 ? max_I : 1.0;
	
	return;
}



void prepare_equations(double A[], double I[], const int N, struct graph *graph_ptr) {
	
	struct node *node_ptr;
	struct edge *edge_ptr;

	int ground_id = graph_ptr->ground->id;
	int source_id = graph_ptr->source->id;
	int id;
	double R;
	
	// for every node check currents going out and comming in
	for(int i=0; i<graph_ptr->node_count; i++) {
		
		// careful with ground and source nodes
		// bit different rules apply to them
		if (i == ground_id || i == source_id)
			continue;
		
		
		// check every edge and make changes to matrix A
		node_ptr = graph_ptr->nodes[i];
		for(uint k=0; k<node_ptr->edges.size(); k++) {
			edge_ptr = node_ptr->edges[k];
			R = edge_ptr->R;
			
			id = edge_ptr->end->id;
			if (id == source_id)
				I[i] -= node_ptr->direction[k] * edge_ptr->end->U / R;
			else if (id != ground_id)
				A[id + N*i] += node_ptr->direction[k] / R;
			
			id = edge_ptr->start->id;
			if (id == source_id)
				I[i] += node_ptr->direction[k] * edge_ptr->start->U / R;
			else if (id != ground_id)
				A[id + N*i] -= node_ptr->direction[k] / R;
		}
	}
	
	// now add to the matrix info about voltage, and source node
	A[graph_ptr->source->id + N*graph_ptr->source->id] = 1;
	I[graph_ptr->source->id] = graph_ptr->source->U;
	
	A[graph_ptr->ground->id + N*graph_ptr->ground->id] = 1;
	I[graph_ptr->ground->id] = 0.0;	//just in case	
	
	return;
}



// 
// double x - value between 0 and 1
// 
// color scaling algorithm derived form site: http://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node76.html
// with changes to scale direction
void calc_RGB_colors(double x, uint8_t *R, uint8_t *G, uint8_t *B) {
	
	x = fmin(1.0, fmax(0.0, x));
	
	// color scale: 0 --> Blue -> Green -> Red --> 1.0
	*R = (uint8_t)fmax(0.0, -510 * fabs(x - 1.0) + 255);
	*G = (uint8_t)fmax(0.0, -510 * fabs(x - 0.5) + 255);
	*B = (uint8_t)fmax(0.0, -510 * fabs(x) + 255);
	
	return;
}




void update_graph(struct graph *graph_ptr, char *file_name) {	
	
	FILE *graph_gv;
	graph_gv = fopen("graph.gv", "w+"); 
	
	
#ifdef DEBUG
	fprintf(graph_gv, "strict digraph circuit {\n"
				      "  node[style=filled, color=pink]\n");
	
	for(int i=0; i<graph_ptr->node_count; i++)
		fprintf(graph_gv, "  %d [label=\"%d: %1.1fV\"]\n",
										graph_ptr->nodes[i]->id,
										graph_ptr->nodes[i]->id,
										graph_ptr->nodes[i]->U );
	
	for(int i=0; i<graph_ptr->edge_count; i++)
		fprintf(graph_gv, "  %d -> %d [color=\"%.1f,1.0,1.0\" label=\"%d, R=%1.1f\\nI=%1.1f\"]\n", 
										graph_ptr->edges[i]->start->id,
										graph_ptr->edges[i]->end->id,
										fabs(graph_ptr->edges[i]->I/graph_ptr->max_I),
										graph_ptr->edges[i]->id,
										graph_ptr->edges[i]->R,
										graph_ptr->edges[i]->I );
										
	fprintf(graph_gv, "  %d [color=\"gray\", label=\"G\"]\n"
					  "  %d [color=\"red\", label=\"S\"]\n"
					  "  %d -> %d [color=\"1.0,1.0,1.0\"]\n"
					  "  %d -> %d [color=\"1.0,1.0,1.0\"]\n"
					  "}",
						graph_ptr->node_count,
						graph_ptr->node_count + 1,
						graph_ptr->node_count,     graph_ptr->ground->id,
						graph_ptr->node_count + 1, graph_ptr->source->id );
	
#else
	fprintf(graph_gv, "strict digraph circuit {\n"
				      "  node[style=filled, color=pink]\n");
	
	uint8_t R, G, B;
	for(int i=0; i<graph_ptr->edge_count; i++) {
		calc_RGB_colors(fabs(graph_ptr->edges[i]->I/graph_ptr->max_I), &R, &G, &B);
		fprintf(graph_gv, "  %d -> %d [label=\"I=%1.2fA\" penwidth=%d color=\"#%02hx%02hx%02hx\"]\n", 
										graph_ptr->edges[i]->start->id,
										graph_ptr->edges[i]->end->id,
										graph_ptr->edges[i]->I,
										(int)(5*fabs(graph_ptr->edges[i]->I/graph_ptr->max_I) + 1),
										R, G, B);
	}
	
	for(int i=0; i<graph_ptr->node_count; i++) {
		calc_RGB_colors(graph_ptr->nodes[i]->U / graph_ptr->source->U, &R, &G, &B);
		fprintf(graph_gv, "  %d [label=\"%d: %1.1fV\" color=\"#%02x%02x%02x\"]\n",
										graph_ptr->nodes[i]->id,
										graph_ptr->nodes[i]->id,
										graph_ptr->nodes[i]->U,
										R, G, B );
	}
	fprintf(graph_gv, "  %d [color=\"gray\", label=\"G\"]\n"
					  "  %d [color=\"red\", label=\"S\"]\n"
					  "  %d -> %d [color=\"1.0,1.0,1.0\"]\n"
					  "  %d -> %d [color=\"1.0,1.0,1.0\"]\n"
					  "}",
						graph_ptr->node_count,
						graph_ptr->node_count + 1,
						graph_ptr->node_count,     graph_ptr->ground->id,
						graph_ptr->node_count + 1, graph_ptr->source->id );
#endif
	
	
	
    fclose(graph_gv);
    
    
    char *command;
    command = (char*)malloc((30 + strlen(file_name)) * sizeof(*command));
    sprintf(command, "sfdp -T png graph.gv > %s.png", file_name);
	system(command);
	free(command);
	
	return;
}




void ini_graph(int node_count, int edge_count, struct graph *graph_ptr) {
	
	graph_ptr->max_I = 1.0;	// put something that is not 0.0
	
	graph_ptr->nodes = new struct node*[node_count];
	graph_ptr->node_count = node_count;
	for(int i=0; i<node_count; i++) {
		graph_ptr->nodes[i] = new struct node;
		graph_ptr->nodes[i]->id = i;
		graph_ptr->nodes[i]->U = 0.0;
	}
	
	
	graph_ptr->edges = new struct edge*[edge_count];
	graph_ptr->edge_count = edge_count;
	for(int i=0; i<edge_count; i++) {
		graph_ptr->edges[i] = new struct edge;
		graph_ptr->edges[i]->id = i;
		graph_ptr->edges[i]->I = 0.0;
	}
	
	return;
}


void load_graph(const char *file, struct graph* graph_ptr) {
	
	//file format:
	//node_count edge_count
	//ground_node source_node voltage
	//node node resistance
	//node node resistance
	//itd...  
	std::ifstream input(file);
	
	
	// allocate space for nodes, edges and initialize some values
	int node_count, edge_count;
	input >> node_count >> edge_count;
	ini_graph(node_count, edge_count, graph_ptr);
	
	// setup ground, source nodes and voltage across them
	int gn, sn;
	double voltage;
	input >> gn >> sn >> voltage;
	graph_ptr->ground = graph_ptr->nodes[gn];
	graph_ptr->source = graph_ptr->nodes[sn];
	graph_ptr->source->U = voltage;
	
	// now load all edges
	int u, v;
	double r;
	for(int i=0; i<edge_count; i++) {
		input >> u >> v >> r;
		
		printf("Adding edge %d -> %d - %.3f ohm\n", u, v, r); 
		
		graph_ptr->edges[i]->R = r;
		graph_ptr->edges[i]->start = graph_ptr->nodes[u];
		graph_ptr->edges[i]->end = graph_ptr->nodes[v];
		
		
		graph_ptr->nodes[u]->edges.push_back(graph_ptr->edges[i]);	//for node u this edge is coming out
		graph_ptr->nodes[u]->direction.push_back(-1);
		
		graph_ptr->nodes[v]->edges.push_back(graph_ptr->edges[i]);	//for node v thie edge is coming in
		graph_ptr->nodes[v]->direction.push_back(1);
	}
	
	input.close();
	return;
}




