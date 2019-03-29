#include <iostream>
#include <utility>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>


struct node{
    int id;
    double U;
    // fields bellow could be integrated into std:pair<>, but this would require changes in code (TODO)
    std::vector<int> direction;		// direction = 1 if edge is coming TO node  OR  -direction = -1 if edge is coming OUT of the node
    std::vector<struct edge*> edges;
};

struct edge{
	int id;	// edge index
	double R;
	double I;
	struct node *start;
	struct node *end;
};


struct graph{
	int node_count;
	int edge_count;
	double max_I;
	struct node *ground;
	struct node *source;
    struct node **nodes;
    struct edge **edges;
};


void ini_graph(int node_count, int edge_count, struct graph *graph_ptr);
void load_graph(const char *file, struct graph* graph_ptr);
void prepare_equations(double A[], double I[], const int N, struct graph *graph_ptr);



void update_graph(struct graph *graph_ptr, char *file_name);
void update_currents(int U[], double I[], struct graph *graph_ptr);

