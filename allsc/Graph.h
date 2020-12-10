#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "LinearHeap.h"

using namespace std;

struct Element{
	Element *next;
	int value;
};

struct Edge {
	Edge *pre, *next;
	Edge *duplicate;
	int adj_id;
};

struct Node {
	Element *head, *tail;
	Edge *first, *last;
};

class Graph {
private:
	string dir;
	int K;
	ui n, m, pm; //pm (previous m) used for progressively reducing the memory useage

	Element *elements; // buffer to allocate all elements
	Edge *edges; // edges of a graph
	Node *nodes; // nodes of a graph

	Node *pnodes; // nodes for super graph
	Edge *pedges; // edges for super graph

	int *computed; //used for marking visited nodes
	int *Q; // queue

	LinearHeap *heap;
	char *inL; // used in decomposition

	int *deg; // used for k-core optimization

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph();
	void compute_allSC() ;

	//void max_spanning_tree() ;
	//void optimization_tree(int output) ;

private:
	void reduce_memory(int nid, int *ids) ;
	void remove_degree_one(FILE *fout) ;
	int kSC_BU(int _K, int *id1, int n_id1, int *id2, int &n_id2, FILE *fout) ;
	void add_edge(Node &node, Edge *edge) ;
	void kCore_prune(int n_Q, FILE *fout);
	int construct_pgraph(int s) ;
	void remove_inter_edges(const vector<Element *> &cc, FILE *fout) ;

	void decomposition(int s, vector<Element *> &cc, int &max_l) ;
	void merge(int s, int t, LinearHeap *heap) ;
	
	void delete_edge(Node &node, Edge *edge) ;
};

#endif
