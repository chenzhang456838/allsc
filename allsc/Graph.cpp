#include "Utility.h"
#include "Graph.h"

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	K = -1;

	nodes = NULL;
	edges = NULL;
	elements = NULL;

	pnodes = NULL;
	pedges = NULL;

	inL = NULL;
	computed = NULL;
	Q = NULL;
	deg = NULL;

	heap = NULL;

	n = m = 0;
}

Graph::~Graph() {
	if(nodes != NULL) delete[] nodes;
	if(edges != NULL) delete[] edges;

	if(pnodes != NULL) delete[] pnodes;
	if(pedges != NULL) delete[] pedges;

	if(inL != NULL) delete[] inL;
	if(computed != NULL) delete[] computed;
	if(Q != NULL) delete[] Q;

	if(elements != NULL) delete[] elements;

	if(deg != NULL) delete[] deg;

	if(heap != NULL) delete heap;
}

void Graph::read_graph() {
	FILE *f = open_file((dir + string("/edges.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	if(nodes == NULL) nodes = new Node[n];
	if(edges == NULL) edges = new Edge[2*m];

	for(ui i = 0;i < n;i ++) nodes[i].first = nodes[i].last = NULL;
	
	if(computed == NULL) computed = new int[n];

	ui n_edges = 0;

	for(ui i = 0;i < n;i ++) {
		int a, d;
		fread(&a, sizeof(int), 1, f);
		fread(&d, sizeof(int), 1, f);
		fread(computed, sizeof(int), d, f);

		for(ui j = 0;j < d;j ++) {
			int b = computed[j];

			edges[n_edges].adj_id = b;
			edges[n_edges].duplicate = &edges[n_edges+1];
			add_edge(nodes[a], &edges[n_edges]);
			++ n_edges;

			edges[n_edges].adj_id = a;
			edges[n_edges].duplicate = &edges[n_edges-1];
			add_edge(nodes[b], &edges[n_edges]);
			++ n_edges;
		}
	}

	fclose(f);
}

void Graph::compute_allSC() {
	if(elements == NULL) elements = new Element[n];
	if(deg == NULL) deg = new int[n];
	if(computed == NULL) computed = new int[n];
	if(Q == NULL) Q = new int[n];

	for(ui i = 0;i < n;i ++) {
		elements[i].value = i;
		computed[i] = 0;
	}

	//write number of nodes and edges (n, m);
	string out_name = dir + "/cg-rsorted.bin";
	FILE *fout = open_file(out_name.c_str(), "wb");
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(int), 1, fout);

	pm = m;

	remove_degree_one(fout);
	printf("number of edges after removing degree-one edges: %u\n", m);

	int *id1 = new int[n], *id2 = new int[n];
	int n_id1, n_id2;

	for(int i = 0;i < n;i ++) id1[i] = i;
	n_id1 = n;

	if(pnodes == NULL) pnodes = new Node[n];
	if(pedges == NULL) pedges = new Edge[2*m];

	//printf("Start computing allsc\n");

	// compute k-edge connected components
	for(int k = 2;;k ++) {
		//printf("%d nodes(%d) edges(%d)\n", k, n_id1, m);
		if(!kSC_BU(k, id1, n_id1, id2, n_id2, fout)) {
			printf("Max SC value: %d\n", k-1);
			break;
		}

		if(pm > 1000000&&m <= pm/2||pm <= 1000000&&m <= pm/8) {
			reduce_memory(n_id2, id2);
			pm = m;
        }

		n_id1 = n_id2;
		int *tmp = id1; id1 = id2; id2 = tmp;
	}
	fclose(fout);

	delete[] id1;
	delete[] id2;

	delete[] nodes; nodes = NULL;
	delete[] edges; edges = NULL;
	delete[] pnodes; pnodes = NULL;
	delete[] pedges; pedges = NULL;
	delete[] elements; elements = NULL;
	delete[] deg; deg = NULL;
	delete[] computed; computed = NULL;

	printf("Starting reverse the cg-rorder.bin file.\n");

	//sort edges in the connectivity graph in decreasing weight order
	FILE *fin = open_file(out_name.c_str(), "rb");
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(int), 1, fin);
	int *buf = new int[3*m];
	fread(buf, sizeof(int), 3*m, fin);

	ui j = 3*m-3;
	for(ui i = 0;i < j;i += 3, j -= 3) {
		swap(buf[i], buf[j]);
		swap(buf[i+1], buf[j+1]);
		swap(buf[i+2], buf[j+2]);
	}

	fclose(fin);

	out_name = dir + "/cg-sorted.bin";
	fout = open_file(out_name.c_str(), "wb");
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(int), 1, fout);
	fwrite(buf, sizeof(int), 3*m, fout);
	fclose(fout);
}

void Graph::reduce_memory(int nid, int *ids) {
	ui edge_c = 0;
	
	for(ui i = 0;i < nid;i ++) {
        for(Edge *e = nodes[ids[i]].first;e != NULL;e = e->next) {
			if(e->adj_id < ids[i]) continue;

			pedges[edge_c].adj_id = e->adj_id;
			pedges[edge_c].duplicate = &pedges[edge_c+1];
			add_edge(pnodes[ids[i]], &pedges[edge_c]);
			++ edge_c;

			pedges[edge_c].adj_id = ids[i];
			pedges[edge_c].duplicate = &pedges[edge_c-1];
			add_edge(pnodes[e->adj_id], &pedges[edge_c]);
			++ edge_c;
		}
        nodes[ids[i]].first = NULL;
	}
    
	delete[] edges;
	edges = new Edge[2*m];

	edge_c = 0;
	for(ui i = 0;i < nid;i ++) {
        for(Edge *e = pnodes[ids[i]].first;e != NULL;e = e->next) {
			if(e->adj_id < ids[i]) continue;

			edges[edge_c].adj_id = e->adj_id;
			edges[edge_c].duplicate = &edges[edge_c+1];
			add_edge(nodes[ids[i]], &edges[edge_c]);
			++ edge_c;

			edges[edge_c].adj_id = ids[i];
			edges[edge_c].duplicate = &edges[edge_c-1];
			add_edge(nodes[e->adj_id], &edges[edge_c]);
			++ edge_c;
		}
        pnodes[ids[i]].first = NULL;
	}

	delete[] pedges;
	pedges = new Edge[2*m];
}

void Graph::remove_degree_one(FILE *fout) {
	ui n_Q = 0;

	if(Q == NULL) Q = new int[n];
	
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		for(Edge *edge = nodes[i].first;edge != NULL;edge = edge->next) ++ cnt;

		if(cnt <= 1) Q[n_Q ++] = i;

		deg[i] = cnt;
	}

	for(int i = 0;i < n_Q;i ++) {
		int s = Q[i];
		computed[s] = 1;

		for(Edge *edge = nodes[s].first;edge != NULL;edge = edge->next) {
			int t = edge->adj_id;

			delete_edge(nodes[t], edge->duplicate);

			if((-- deg[t]) == 1) Q[n_Q ++] = t;

			int tt[3];
			tt[0] = s; tt[1] = t; tt[2] = 1;
			if(tt[0] > tt[1]) swap(tt[0], tt[1]);

			fwrite(tt, sizeof(int), 3, fout);
			-- m;
		}

		nodes[s].last = nodes[s].first = NULL;
	}
}

int Graph::kSC_BU(int _K, int *id1, int n_id1, int *id2, int &n_id2, FILE *fout) {
	K = _K;

	for(ui i = 0;i < n_id1;i ++) computed[id1[i]] = 0;

	ui n_Q = 0;
	for(ui i = 0;i < n_id1;i ++) {
		int cnt = 0;
		for(Edge *edge = nodes[id1[i]].first;edge != NULL;edge = edge->next) ++ cnt;

		if(cnt < K) Q[n_Q ++] = id1[i];

		deg[id1[i]] = cnt;
	}

	kCore_prune(n_Q, fout);

	int max_l = 0, non_trivial = 0;

	n_id2 = 0;
	for(int i = 0;i < n_id1;i ++) {
		// if(i%10000 == 0) printf(".");

		if(computed[id1[i]]) continue;
        
		if(construct_pgraph(id1[i]) > 1) non_trivial = 1;

		vector<Element *> cc;
		decomposition(id1[i], cc, max_l);

		if(cc.size() == 1) {
			for(Element *e = cc[0];e != NULL;e = e->next) {
				computed[e->value] = 1;
				id2[n_id2 ++] = e->value;
			}
		}
		else remove_inter_edges(cc, fout);

		-- i;
	}

	return non_trivial;
}

int Graph::construct_pgraph(int s) {
	ui pedge_c = 0;
	int q_c = 1;

	computed[s] = 1; Q[0] = s;

	for(int i = 0;i < q_c;i ++) {
		s = Q[i];

		pnodes[s].head = pnodes[s].tail = &elements[s];
		elements[s].next = NULL;

		for(Edge *edge = nodes[s].first;edge != NULL;edge = edge->next) {
			if(!computed[edge->adj_id]) {
				computed[edge->adj_id] = 1;
				Q[q_c ++] = edge->adj_id;
			}

			if(edge->adj_id > s) {
				int a = s, b = edge->adj_id;

				pedges[pedge_c].adj_id = b;
				pedges[pedge_c].duplicate = &pedges[pedge_c+1];
				add_edge(pnodes[a], &pedges[pedge_c]);
				++ pedge_c;

				pedges[pedge_c].adj_id = a;
				pedges[pedge_c].duplicate = &pedges[pedge_c-1];
				add_edge(pnodes[b], &pedges[pedge_c]);
				++ pedge_c;
			}
		}
	}

	for(ui i = 0;i < q_c;i ++) computed[Q[i]] = 0;

	return q_c;
}

void Graph::kCore_prune(int n_Q, FILE *fout) {
	for(int i = 0;i < n_Q;i ++) {
		int s = Q[i];
		computed[s] = 1;

		for(Edge *edge = nodes[s].first;edge != NULL;edge = edge->next) {
			int t = edge->adj_id;

			delete_edge(nodes[t], edge->duplicate);

			if((-- deg[t]) == K-1) Q[n_Q ++] = t;

			int tt[3];
			tt[0] = s; tt[1] = t; tt[2] = K-1;
			if(tt[0] > tt[1]) swap(tt[0], tt[1]);

			fwrite(tt, sizeof(int), 3, fout);
			-- m;
		}

		nodes[s].last = nodes[s].first = NULL;
	}
}

void Graph::merge(int s, int t, LinearHeap *heap) {
	pnodes[s].tail->next = pnodes[t].head;
	pnodes[s].tail = pnodes[t].tail;

	Edge *e = pnodes[t].first;
	Edge *tmp;

	while(e != NULL) {
		tmp = e->next;
        
		if(e->adj_id == s) {
			if(heap != NULL) heap->set_key(s, heap->get_key(s) - 1);
			delete_edge(pnodes[e->adj_id], e->duplicate);
		}
		else {
			e->duplicate->adj_id = s;
			add_edge(pnodes[s], e);
		}

		e = tmp;
	}

	pnodes[t].first = NULL;
}

void Graph::decomposition(int ss, vector<Element *> &cc, int &max_l) {
	if(heap == NULL) heap = new LinearHeap(n);

	if(inL == NULL) {
		inL = new char[n];
		memset(inL, 0, sizeof(char)*n);
	}

	cc.clear();

	int cnt = 0;

	while(pnodes[ss].first != NULL) {
		int s = ss;

		++ cnt;

		heap->insert(s, 0);

		int n_Q = 0;
		int key;

		while(1) {
			if(!heap->extract_max(s, key)) break;

			inL[s] = 1;
			Q[n_Q ++] = s;

			int new_qc = n_Q;
			for(int i = n_Q - 1;i < new_qc;i ++) {
				int u = Q[i];
                for(Edge *e = pnodes[u].first;e != NULL;e = e->next) {
                    if(!inL[e->adj_id]) {
                        int new_key = heap->get_key(e->adj_id);
                        if(new_key < K) {
                            if(new_key > 0) heap->remove(e->adj_id);

                            new_key += 1;
                            if(new_key >= K) {
                                heap->set_key(e->adj_id, new_key);
                                Q[new_qc ++] = e->adj_id;
                            }
                            else heap->insert(e->adj_id, new_key);
                        }
                        else heap->set_key(e->adj_id, new_key + 1);
                    }
                }
                
				if(u == s) continue;

				heap->set_key(s, heap->get_key(s) + heap->get_key(u));
				heap->set_key(u, 0);
				inL[u] = 0;
				merge(s, u, heap);
			}
		}
        
        -- n_Q;
		while(n_Q > 0&&heap->get_key(Q[n_Q]) < K) {
			int t = Q[n_Q]; -- n_Q;
			cc.push_back(pnodes[t].head);

			heap->set_key(t, 0);
			inL[t] = 0;

			for(Edge *e = pnodes[t].first;e != NULL;e = e->next) delete_edge(pnodes[e->adj_id], e->duplicate);

			pnodes[t].first = NULL;
		}

		for(int i = 0;i <= n_Q;i ++) {
			heap->set_key(Q[i], 0);
			inL[Q[i]] = 0;
		}
	}

	if(cnt > max_l) max_l = cnt;

	cc.push_back(pnodes[ss].head);
}

void Graph::remove_inter_edges(const vector<Element *> &cc, FILE *fout) {
	for(int j = 0;j < (int)cc.size();j ++) for(Element *e = cc[j];e != NULL;e = e->next) computed[e->value] = j+1;
	
	int q_c = 0;

	for(int j = 0;j < (int)cc.size();j ++) for(Element *e = cc[j];e != NULL;e = e->next) {
		int s = e->value;
		Edge *list = nodes[s].first;
		nodes[s].first = nodes[s].last = NULL;
		int cnt = 0;

		while(list != NULL) {
			Edge *tmp = list->next;
			if(computed[list->adj_id] == computed[s]) {
				if(nodes[s].first == NULL) {
					nodes[s].first = nodes[s].last = list;
					list->pre = NULL;
				}
				else {
					nodes[s].last->next = list;
					list->pre = nodes[s].last;
					nodes[s].last = list;
				}

				++ cnt;
			}
			else if(s < list->adj_id) {
				int tt[3];
				tt[0] = s; tt[1] = list->adj_id; tt[2] = K-1;
				fwrite(tt, sizeof(int), 3, fout);
				-- m;
			}

			list = tmp;
		}

		if(nodes[s].last != NULL) nodes[s].last->next = NULL;

		deg[s] = cnt;
		if(cnt < K) Q[q_c ++] = s;
	}

	for(int j = 0;j < (int)cc.size();j ++) for(Element *e = cc[j];e != NULL;e = e->next) computed[e->value] = 0;

	kCore_prune(q_c, fout);
}

void Graph::add_edge(Node &node, Edge *edge) {
	edge->next = NULL;

	if(node.first == NULL) {
		node.first = node.last = edge;
		edge->pre = NULL;
	}
	else {
		node.last->next = edge;
		edge->pre = node.last;
		node.last = edge;
	}
}

void Graph::delete_edge(Node &node, Edge *edge) {
	if(edge->pre == NULL) {
		node.first = edge->next;
		if(edge->next != NULL) edge->next->pre = NULL;
	}
	else {
		if(edge == node.last) node.last = edge->pre;

		Edge *tmp = edge->pre;
		tmp->next = edge->next;

		if(edge->next != NULL) edge->next->pre = tmp;
	}
}
