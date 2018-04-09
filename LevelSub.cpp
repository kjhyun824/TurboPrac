#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <stack>
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#include "LevelSub.h"
#include "StringUtility.h"
#include "runtimecounter.h"

LevelSub::LevelSub(const std::string DataFileName, const std::string QueryFileName, const int topk, const int debug) :
	DataGraphFileName(DataFileName), QueryGraphFileName(QueryFileName), k(topk)
{
	if(debug > 0)
		debugflag = true;
	LoadGraph();
	LoadQuery();
}

LevelSub::~LevelSub()
{

}

/* KJH - Rewrite Query graph to NEC Tree */
void LevelSub::rewriteToNECTree() {
	Graph q = q_list[0];
	dirGraph NEC_tree;
	NEC_vertex_map_t NEC_Data;
	std::map<vertex_t,bool> query_visited; // visited_query_nodes

	vertex_iterator_t v_bg, v_end;
	tie(v_bg,v_end) = vertices(q);

	/* Set the root node */
	NEC_Data[0].push_back(*v_bg);
	dir_vertex_t root = add_vertex(*v_bg,NEC_tree);
	NEC_vertices.push_back(root);
	NEC_label[root] = label_q_list[0][*v_bg];

	query_visited[*v_bg++] = true; // mark root visited
	for(; v_bg != v_end; v_bg++) {
		query_visited[*v_bg] = false; 
	}

	std::vector<dir_vertex_t> v_curr, v_next, tmp;
	v_next.push_back(root);


	dir_vertex_t added_v;
	while(!v_next.empty()) {
		/*
		while(!v_next.empty()) { // v_curr = v_next , v_next = empty
			v_curr.push_back(v_next.front());
			v_next.erase(v_next.begin());
		}
		*/

		/* v_curr = v_next , v_next = empty */
		tmp = v_curr;
		v_curr = v_next;
		v_next = tmp;
		v_next.clear();

		for(std::vector<dir_vertex_t>::iterator nec_vid = v_curr.begin(); nec_vid != v_curr.end(); ++nec_vid) { // iterations of NEC vertices
			std::set<vertex_t> neighbors;
			std::set<int> label_set;
			std::map<int,std::vector<vertex_t> > group_by_label;

			adjacency_iterator_t adj_t, adj_end;
			
			for(std::vector<vertex_t>::iterator nec_it = NEC_Data[*nec_vid].begin(); nec_it != NEC_Data[*nec_vid].end(); ++nec_it) {
				for(tie(adj_t,adj_end) = adjacent_vertices(*nec_it,q); adj_t != adj_end; ++adj_t) {
					if(!query_visited[*adj_t]) {
						neighbors.insert(*adj_t);
						query_visited[*adj_t] = true;
					}
				}
			}

			if(!neighbors.empty()) {
				for(std::set<vertex_t>::iterator v_it = neighbors.begin(); v_it != neighbors.end(); ++v_it) {
					label_set.insert(label_q_list[0][*v_it]);
					group_by_label[label_q_list[0][*v_it]].push_back(*v_it);
				}

				for(std::set<int>::iterator label_it = label_set.begin(); label_it != label_set.end(); ++label_it) {
					std::vector<std::vector<vertex_t> > NECVertices = FindNEC(group_by_label[*label_it]);
					for(std::vector<std::vector<vertex_t> >::iterator nec_v_it = NECVertices.begin(); nec_v_it != NECVertices.end(); ++nec_v_it) {
						int nec_idx = num_vertices(NEC_tree);

						for(std::vector<vertex_t>::iterator vertex = (*nec_v_it).begin(); vertex != (*nec_v_it).end(); ++vertex)
							NEC_Data[nec_idx].push_back(*vertex);

						added_v = add_vertex(nec_idx,NEC_tree);
						NEC_vertices.push_back(added_v);
						add_edge(*nec_vid,added_v,NEC_tree);
						NEC_label[added_v] = *label_it;
						v_next.push_back(nec_idx);
						parent_map[added_v] = *nec_vid;
					}
				}
			}
		}
	}

	NEC_data_list.push_back(NEC_Data);
	NEC_tree_list.push_back(NEC_tree);

	std::pair<dir_vertex_iterator_t, dir_vertex_iterator_t> vp;
	for (vp = vertices(NEC_tree); vp.first != vp.second; ++vp.first) {
		printf("[DEBUG] NEC_Data | vid : %d(label:%d), vertices : ",(int) *vp.first,NEC_label[*vp.first]);
		for(int i = 0; i < NEC_Data[*vp.first].size(); i++) 
			printf("%d ",(int) NEC_Data[*vp.first][i]);
		printf("\n");
	}

	dir_edge_iterator_t ei, ei_end;
	printf("[DEBUG] edges(g) = ");
	for (tie(ei, ei_end) = edges(NEC_tree); ei != ei_end; ++ei)
		printf("(%d,%d) ",(int) source(*ei, NEC_tree),(int) target(*ei, NEC_tree));
	printf("\n");
}
/* End of rewriteToNECTree */

/* KJH - FindNEC Function */
std::vector<std::vector<vertex_t> > LevelSub::FindNEC(std::vector<vertex_t> v) {
	Graph q = q_list[0];
	std::vector<std::vector<vertex_t> > new_NEC;

	if(v.size() == 1) {
		new_NEC.push_back(v);
		return new_NEC;
	}

	std::set<int> adj_set[v.size()];

	bool is_alive[v.size()];
	for(int i = 0; i < v.size(); i++) is_alive[i] = true;

	for(int i = 0; i < v.size(); i++) {
		adjacency_iterator_t adj_t, adj_end;
		for(tie(adj_t, adj_end) = adjacent_vertices(v[i], q); adj_t != adj_end; ++adj_t) 
			adj_set[i].insert(*adj_t);
	}

	for(int i = 0; i < v.size(); i++) {
		bool inserted = false;
		std::vector<vertex_t> temp;
		if(!is_alive[i]) continue;

		for(int j = i+1; j < v.size(); j++) {
			if(!is_alive[j]) {
				break;
			}

			if(adj_set[i] == adj_set[j]) {
				if(!inserted) {
					temp.push_back(v[i]);
					new_NEC.push_back(temp);
					inserted = true;
				}
				new_NEC[new_NEC.size()-1].push_back(v[j]);
				is_alive[j] = false;
			}
		}

		if(!inserted) { //If the class has a unique vertex
			temp.push_back(v[i]);
			new_NEC.push_back(temp);
		}
	}

	return new_NEC;
}

bool LevelSub::NLF(vertex_t u, vertex_t v) {
	Graph q = q_list[0];
	std::map<int, int> query_label_count, data_label_count;
	adjacency_iterator_t adj_t, adj_end;

	int curr_label;

	for(tie(adj_t,adj_end) = adjacent_vertices(u,q); adj_t != adj_end; ++adj_t) {
		curr_label = label_q_list[0][*adj_t];
		if(query_label_count.find(curr_label) == query_label_count.end()) query_label_count[curr_label] = 1;
		else query_label_count[curr_label] += 1;
	}

	for(tie(adj_t,adj_end) = adjacent_vertices(v,g); adj_t != adj_end; ++adj_t) {
		curr_label = label_v[*adj_t];
		if(data_label_count.find(curr_label) == data_label_count.end()) data_label_count[curr_label] = 1;
		else data_label_count[curr_label] += 1;
	}

	for(std::map<int,int>::iterator it = query_label_count.begin(); it != query_label_count.end(); ++it) {
		if(data_label_count.find(it->first) != data_label_count.end() && data_label_count[it->first] >= query_label_count[it->first]) continue;
		else return false;
	}

	return true;
}


dir_vertex_t nec_tree_child(dir_vertex_t u, int idx) {
}

bool LevelSub::exploreCR(dir_vertex_t u, std::vector<vertex_t> V, vertex_t v_parent) {
	Graph q = q_list[0];
	dirGraph NEC_tree = NEC_tree_list[0];
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	std::vector<dir_vertex_t> children;
	std::vector<int> child_label_list;

	for(int v_idx = 0; v_idx < V.size(); v_idx++) {
		if(data_visited[V[v_idx]]) continue;

		if(out_degree(NEC_Data[u][0],q) > out_degree(V[v_idx],g) || !NLF(NEC_Data[u][0],V[v_idx])) {
			continue;
		}

		data_visited[V[v_idx]] = true;
		bool matched = true;

		std::vector<dir_vertex_t> sorted_children;
		std::vector<int> adj_size;
		std::map<dir_vertex_t,std::vector<vertex_t> > v_map_adj_list;
		std::map<int,std::vector<dir_vertex_t> > size_map_vertex;

		dir_adjacency_iterator_t dir_adj_t,dir_adj_end;
		for(tie(dir_adj_t,dir_adj_end) = adjacent_vertices(u,NEC_tree); dir_adj_t != dir_adj_end; ++dir_adj_t) {
			int adj_count = 0;
			adjacency_iterator_t adj_t,adj_end;
			for (tie(adj_t, adj_end) = adjacent_vertices(V[v_idx], g); adj_t != adj_end; ++adj_t) {
				if(label_v[*adj_t] == NEC_label[*dir_adj_t]) {
					v_map_adj_list[*dir_adj_t].push_back(*adj_t);
					adj_count++;
				}
			}
			adj_size.push_back(adj_count);
			size_map_vertex[adj_count].push_back(*dir_adj_t);
		}

		std::sort (adj_size.begin(),adj_size.end());

		for(int i = 0; i < adj_size.size(); i++) {
			while(!size_map_vertex[adj_size[i]].empty()) {
				sorted_children.push_back(size_map_vertex[adj_size[i]].front());
				size_map_vertex[adj_size[i]].erase(size_map_vertex[adj_size[i]].begin());
			}
		}
		/* Sorting Done */

		for(int i = 0; i < sorted_children.size(); i++) {
			if(!exploreCR(sorted_children[i],v_map_adj_list[sorted_children[i]],V[v_idx])) {
				for(int j = 0; j < sorted_children.size() && j != i; j++) 
					clearCR(sorted_children[j],V[v_idx]);
				matched = false;
				break;
			}
		}

		data_visited[V[v_idx]] = false;
		if(!matched) continue;

		cand_sub[std::make_pair(u,v_parent)].push_back(V[v_idx]);
		cand_key_pair[u].insert(v_parent);
	}

	if(cand_sub[std::make_pair(u,v_parent)].size() < NEC_Data[u].size()) {
		clearCR(u,v_parent);
		return false;
	}

	return true;
}

void LevelSub::clearCR(dir_vertex_t u, vertex_t v) {
	cand_sub[std::make_pair(u,v)].clear();
	cand_sub.erase(std::make_pair(u,v));
	cand_key_pair[u].erase(v);
	if(cand_key_pair[u].empty()) cand_key_pair.erase(u);
}

int LevelSub::cand_sub_size(dir_vertex_t u) {
	int count = 0;
	std::set<vertex_t> pair_keys = cand_key_pair[u];
	for(std::set<vertex_t>::iterator it = pair_keys.begin(); it != pair_keys.end(); ++it) {
		count += cand_sub[std::make_pair(u,*it)].size();
	}

	return count;
}

double LevelSub::non_tree_edge_num(dir_vertex_t u) {
	int count = 0;

	Graph q = q_list[0];
	dirGraph NEC_tree = NEC_tree_list[0];
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	dir_vertex_t curr_parent = parent_map[u];
	vertex_t curr_v = NEC_Data[u][0], parent_v = NEC_Data[curr_parent][0];
	
	adjacency_iterator_t adj_t, adj_end;
	for(tie(adj_t,adj_end) = adjacent_vertices(curr_v,q); adj_t != adj_end; ++adj_t) {
		//parent?
		bool is_parent = false;
		for(std::vector<vertex_t>::iterator it = NEC_Data[curr_parent].begin(); it != NEC_Data[curr_parent].end(); ++it) 
			if(*adj_t == *it) {
				is_parent = true;
				break;
			}
		
		if(is_parent) continue;
		//child?
		bool is_child = false;
		dir_adjacency_iterator_t dir_adj_t,dir_adj_end;
		for(tie(dir_adj_t,dir_adj_end) = adjacent_vertices(u,NEC_tree); dir_adj_t != dir_adj_end; ++dir_adj_t) {
			for(std::vector<vertex_t>::iterator it = NEC_Data[*dir_adj_t].begin(); it != NEC_Data[*dir_adj_t].end(); ++it)
				if(*adj_t == *it) {
					is_child = true;
					break;
				}
			if(is_child) break;
		}

		if(is_child) continue;

		count++;
	}

	u = curr_parent;
	curr_parent = parent_map[curr_parent];
	for(tie(adj_t,adj_end) = adjacent_vertices(parent_v,q); adj_t != adj_end; ++adj_t) {
		//parent?
		bool is_parent = false;
		for(std::vector<vertex_t>::iterator it = NEC_Data[curr_parent].begin(); it != NEC_Data[curr_parent].end(); ++it) 
			if(*adj_t == *it) {
				is_parent = true;
				break;
			}
		
		if(is_parent) continue;
		//child?
		bool is_child = false;
		dir_adjacency_iterator_t dir_adj_t,dir_adj_end;
		for(tie(dir_adj_t,dir_adj_end) = adjacent_vertices(u,NEC_tree); dir_adj_t != dir_adj_end; ++dir_adj_t) {
			for(std::vector<vertex_t>::iterator it = NEC_Data[*dir_adj_t].begin(); it != NEC_Data[*dir_adj_t].end(); ++it)
				if(*adj_t == *it) {
					is_child = true;
					break;
				}
			if(is_child) break;
		}

		if(is_child) continue;

		count++;
	}

	return ((double) count);
}

std::vector<dir_vertex_t> LevelSub::find_path(dir_vertex_t u) {
	std::vector<dir_vertex_t> path;
	path.push_back(u);

	if(u == 0) return path;

	dir_vertex_t u_next = u;
	while(u_next != 0) {
		u_next = parent_map[u_next];
		path.push_back(u_next);
	}

	std::reverse(path.begin(),path.end());

	return path;
}

void LevelSub::matchingOrder() {
	dirGraph NEC_tree = NEC_tree_list[0];
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	std::vector<dir_vertex_t> leafnodes;
	std::vector<double> leafprice;
	dir_vertex_iterator_t v_bg,v_end;

	for(tie(v_bg,v_end) = vertices(NEC_tree); v_bg != v_end; ++v_bg) {
		if(out_degree(*v_bg,NEC_tree) == 0) {
			leafnodes.push_back(*v_bg);
			double price = (((double) cand_sub_size(*v_bg)) / (non_tree_edge_num(*v_bg) + 1.0));
			leafprice.push_back(price);
			printf("[DEBUG] for %d, price is %.2f\n",(int) *v_bg,price);
		}
	}

	std::vector<dir_vertex_t> matching_order;
	while(!leafnodes.empty()) {
		int min_index;
		double min_price = INT_MAX;
		for(int i = 0; i < leafnodes.size(); i++) 
			if(min_price >= leafprice[i]) {
				min_index = i;
				min_price = leafprice[i];
			}

		std::vector<dir_vertex_t> path_of_min = find_path(leafnodes[min_index]);

		while(!path_of_min.empty()) {
			if(std::find(matching_order.begin(),matching_order.end(),path_of_min.front()) != matching_order.end()) {
				path_of_min.erase(path_of_min.begin());
				continue;
			}

			matching_order.push_back(path_of_min.front());
			path_of_min.erase(path_of_min.begin());
		}

		leafnodes.erase(leafnodes.begin()+min_index);
		leafprice.erase(leafprice.begin()+min_index);
	}

	matching_order_list.push_back(matching_order);
}

std::vector<vertex_t> LevelSub::NextComb(dir_vertex_t u,std::vector<vertex_t> C) {
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	int nec_size = NEC_Data[u].size();
	std::vector<vertex_t> next_comb;

	if(C.size() < nec_size) 
		printf("[DEBUG] Candidate size is smaller than NEC_data size\n");

	if(curr_comb_idx_map.empty() || curr_comb_idx_map.find(u) == curr_comb_idx_map.end()) {
		for(int i = 0; i < nec_size; i++) {
			curr_comb_idx_map[u].push_back(i);
			next_comb.push_back(C[i]);
		}

		curr_comb_map[u] = next_comb;

		return next_comb;
	}

	if(curr_comb_idx_map[u][0] >= C.size() - nec_size) {
		printf("[DEBUG] Final Comb is Done\n");
		curr_comb_idx_map.erase(u);
		return next_comb;
	}

	bool temp_updated = false;
	int comb_iter = curr_comb_idx_map[u].size() - 1;

	while(!temp_updated) {
		curr_comb_idx_map[u][comb_iter]++;
		if(curr_comb_idx_map[u][comb_iter] == C.size() - (nec_size - comb_iter - 1)) {
			comb_iter--;
			continue;
		}

		temp_updated = true;
		comb_iter++;
	}

	for(;comb_iter < nec_size; comb_iter++) {
		curr_comb_idx_map[u][comb_iter] = curr_comb_idx_map[u][comb_iter-1] + 1;
	}

	for(int i = 0; i < curr_comb_idx_map[u].size(); i++) {
		next_comb.push_back(C[curr_comb_idx_map[u][i]]);
	}

	curr_comb_map[u] = next_comb;

	return next_comb;
}

void LevelSub::GenPerm(int idx) {
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	if(idx == NEC_vertices.size()) return;

	if(NEC_Data[idx].size() == 1) {
		GenPerm(idx+1);
	}
	else {
		while(!NextPerm(NEC_vertices[idx]))	{
			GenPerm(idx+1);
		}
	}
}

bool LevelSub::NextPerm(dir_vertex_t u_i) {
	NEC_vertex_map_t NEC_Data = NEC_data_list[0];
	int nec_size = NEC_Data[u_i].size();
	std::vector<vertex_t> curr_comb = curr_comb_map[u_i];
		
	if(curr_perm_idx_map.empty() || curr_perm_idx_map.find(u_i) == curr_perm_idx_map.end()) {
		for(int i = 0; i < nec_size; i++) {
			curr_perm_idx_map[u_i].push_back(i);
			map_table[NEC_Data[u_i][i]] = curr_comb_map[u_i][i];
		}

		return true;
	}

	bool isDone = true;

	for(int i = 0; i < nec_size; i++) {
		if(curr_perm_idx_map[u_i][i] == nec_size-i-1) continue;

		isDone = false;
		break;
	}

	if(isDone) {
		printf("[DEBUG] All permutation is Done\n");
		curr_perm_idx_map[u_i].clear();
		curr_perm_idx_map.erase(u_i);
		return false;
	}

	int temp_arr[nec_size];
	for(int i = 0; i < nec_size; i++) {
		temp_arr[i] = curr_perm_idx_map[u_i].front();
		curr_perm_idx_map[u_i].erase(curr_perm_idx_map[u_i].begin());
	}

	std::next_permutation(temp_arr,temp_arr+nec_size);


	for(int i = 0; i < nec_size; i++) {
		curr_perm_idx_map[u_i].push_back(temp_arr[i]);
		map_table[NEC_Data[u_i][i]] = curr_comb_map[u_i][temp_arr[i]];
	}

	return true;
}

bool LevelSub::IsJoinable(vertex_t v, vertex_t C_i) {
}

void LevelSub::SubgraphSerach() {
}


inline void LevelSub::buildSignature(int u_label, int v, int u)
{
	neighborhoodSignature[v][u_label].push_back(u);
}

void LevelSub::LoadGraph()
{
	std::ifstream dataGraphFile = std::ifstream(DataGraphFileName);
	if (!dataGraphFile.is_open()){
		std::cout << "data graph file doesn't exist" << std::endl;
		exit(1);
	}
	std::string line;
	std::vector<int> integerValues;
	std::getline(dataGraphFile, line);

	while (line.size() != 0 && (*line.begin()) == 't'){
		int number_of_vertices;
		int number_of_edges;
		int number_of_labels;
		std::cout << "Start to read the graph..." << std::endl;
		if (std::getline(dataGraphFile, line) && (*line.begin() == '%'))
		{
			String_Utility::readIntegersFromString(line, integerValues);
			number_of_vertices = integerValues[0];
			number_of_edges = integerValues[1];
			number_of_labels = integerValues[2];
		}

		std::map<int, int> label;
		candidates = std::vector<std::vector<int> >(number_of_labels+1, std::vector<int>());
		neighborhoodSignature.resize(number_of_vertices);

		std::cout << "reading vertices..." << std::endl;
		while (std::getline(dataGraphFile, line) && (*line.begin() == 'v'))
		{
			String_Utility::readIntegersFromString(line, integerValues);
			vertex_t v = add_vertex(integerValues[0], g);
			candidates[integerValues[1]].push_back(v);
			label[v] = integerValues[1];

			/* KJH */
			data_visited[v] = false;
			label_v[v] = integerValues[1];
		}
		std::cout << "reading edges..." << std::endl;
		while ((*line.begin()) == 'e')
		{
			String_Utility::readIntegersFromString(line, integerValues);

			std::pair<vertex_t, vertex_t> newedge(integerValues[0], integerValues[1]);
			add_edge(integerValues[0], integerValues[1], g);
			buildSignature(label[integerValues[1]], integerValues[0], integerValues[1]);
			buildSignature(label[integerValues[0]], integerValues[1], integerValues[0]);

			if (!std::getline(dataGraphFile, line))
				break;
		}

		std::cout << "done! Edges:" << num_edges(g) << " Vertices:" << num_vertices(g) << std::endl;
	}
}

void LevelSub::LoadQuery()
{
	std::vector<int> label_q;
	std::ifstream qGraphFile = std::ifstream(QueryGraphFileName);
	if (!qGraphFile.is_open()){
		std::cout << "data graph file doesn't exist" << std::endl;
		exit(1);
	}
	std::string line;
	std::vector<int> integerValues;
	std::getline(qGraphFile, line);
	std::cout << "Start to read the query graph file..." << std::endl;
	while (line.size() != 0 && (*line.begin()) == 't'){
		int number_of_vertices;
		int number_of_edges;
		int number_of_labels;
		if (std::getline(qGraphFile, line) && (*line.begin() == '%'))
		{
			String_Utility::readIntegersFromString(line, integerValues);
			number_of_vertices = integerValues[0];
			number_of_edges = integerValues[1];
			number_of_labels = integerValues[2];
		}

		Graph q;
		label_q.resize(number_of_vertices);
		while (std::getline(qGraphFile, line) && (*line.begin() == 'v'))
		{
			String_Utility::readIntegersFromString(line, integerValues);
			vertex_t v = add_vertex(integerValues[0], q);
			label_q[v] = integerValues[1];
		}
		while ((*line.begin()) == 'e')
		{
			String_Utility::readIntegersFromString(line, integerValues);
			add_edge(integerValues[0], integerValues[1], q);
			if (!std::getline(qGraphFile, line))
				break;
		}
		q_list.push_back(q);
		label_q_list.push_back(label_q);
	}
}

void LevelSub::randomShuffle()
{
	 // obtain a time-based seed:
	 unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	 //shuffle neighbors
	 for(int i = 0; i < neighborhoodSignature.size(); ++i)
	 {
		 for(auto n_t = neighborhoodSignature[i].begin(); n_t != neighborhoodSignature[i].end(); ++n_t)
		 {
			 shuffle( n_t->second.begin(), n_t->second.end(), std::default_random_engine(seed) );
		 }
	 }

	 //shuffle candidates
	 for(auto cand_t = candidates.begin(); cand_t != candidates.end(); ++cand_t)
	 {
		 shuffle( cand_t->begin(), cand_t->end(), std::default_random_engine(seed) );
	 }
}

void LevelSub::queryOrder(Graph &q, std::vector<int> &label_q, std::vector<int>& sorted)
{
	std::vector<std::pair<int, double> > rank(num_vertices(q), std::make_pair(0, 0.0));
	vertex_iterator_t v_bg, v_end;
	for (tie(v_bg, v_end) = vertices(q); v_bg != v_end; ++v_bg)
	{
		rank[*v_bg] = std::make_pair(*v_bg, static_cast<double>(candidates[label_q[*v_bg]].size()) / static_cast<double>(degree(*v_bg, q)));
	}

	//simply sort all the nodes
	sort(rank.begin(), rank.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b) { return a.second < b.second; });
	for_each(rank.begin(), rank.end(), [&sorted](const std::pair<int, double> &b) { sorted.push_back(b.first); });
}


void LevelSub::subgraphSearch()
{

//	randomShuffle();

	rewriteToNECTree();

	std::vector<vertex_t> temp;
	temp.push_back(0);
	exploreCR(0,temp,-1);

	for(my_it iter = cand_sub.begin(); iter != cand_sub.end(); ++iter) {
		printf("[DEBUG] cand_sub(%d,%d) - ",(int) (*iter).first.first,(int) (*iter).first.second);
		for(int i = 0; i < cand_sub[iter->first].size(); i++)
			printf("%d ",(int) cand_sub[iter->first][i]);
		printf("\n");
	}
	
	matchingOrder();
	std::vector<dir_vertex_t> mo = matching_order_list[0];
	printf("[DEBUG] matching order : ( ");
	for(std::vector<dir_vertex_t>::iterator it = mo.begin(); it != mo.end(); ++it)
		printf("%d ",(int) *it);
	printf(")\n");

	std::vector<vertex_t> temp_C;
	temp_C.push_back(2);
	temp_C.push_back(3);
	temp_C.push_back(4);
	for(int i = 0; i < 3; i++) {
		std::vector<vertex_t> next = NextComb(2,temp_C);

		if(next.empty()) continue;

		printf("[DEBUG] %d - comb : ",i+1);
		for(std::vector<vertex_t>::iterator it = next.begin(); it != next.end(); ++it) {
			printf("%d ",(int) *it);
		}
		printf("\n");

		for(int j = 0; j < 3; j++) {
			if(!NextPerm(2)) continue;
			printf("[DEBUG] Perm : ");
			for(int k = 0; k < next.size(); k++) {
				printf("%d ",(int) curr_comb_map[2][curr_perm_idx_map[2][k]]);
			}
			printf("\n");
		}
	}

	for (int i = 0; i < q_list.size(); ++i)
	{
		matched = d_bitset(num_vertices(g));
		embeddingList.clear();
		nsQ.clear();
		gCT.clear();

		//tmp candidates for every query
		std::vector<std::vector<int> > tmpcand(candidates.begin(), candidates.end());

		Runtimecounter timer; timer.start(); auto gct = std::vector<bool>(num_vertices(q_list[i]), false); gCT = std::vector<std::vector<bool> >(num_vertices(q_list[i]), gct);

		std::map<int, int> labelnumber;
		//ns and gCT for query
		nsQ.resize(num_vertices(q_list[i]));
		vertex_iterator_t q_t, q_end;
		tie(q_t, q_end) = vertices(q_list[i]);
		std::set<vertex_t> query_set(q_t, q_end);
		for (; q_t != q_end; ++q_t)
		{
			std::map<int, std::set<int> > ns;
			adjacency_iterator_t adj_t, adj_end;
			for (tie(adj_t, adj_end) = adjacent_vertices((*q_t), q_list[i]); adj_t != adj_end; ++adj_t)
			{
				int label = label_q_list[i][*adj_t];
				if (ns.find(label) == ns.end())
					ns[label] = std::set<int>();
				ns[label].insert(*adj_t);
				gCT[*q_t][*adj_t] = true;
			}
			nsQ[*q_t] = ns;
			int labelq = label_q_list[i][*q_t];
			if(labelnumber.find(labelq) == labelnumber.end())
				labelnumber[labelq] = 0;
			++labelnumber[labelq];
		}

		std::vector<int> sorted;
		queryOrder(q_list[i], label_q_list[i], sorted);


		std::vector<int> rsorted(sorted.size(), -1);
		for (int j = 0; j < sorted.size(); ++j)
			rsorted[sorted[j]] = j;

		std::vector<int> embedding(num_vertices(q_list[i]), -1);

		//re sort the query
		int overlapSize = 0;
		std::vector<std::pair<int, int> > nSorted;

		chosen = std::vector<bool>(sorted.size(), false);
		chosen[sorted[0]] = true;
		std::vector<bool> noChosenNb(sorted.size(), false);
		//store the iterator of data sets for chosen overlap nodes
		//getPairSorted(overlapQ, sorted, nSorted, q_list[i]);
		lbnum = labelnumber;
		lbRm = std::vector<int>(sorted.size(), 0);
		nbRm = std::vector<int>(sorted.size(), 0);
		getPairSorted3(rsorted, sorted, nSorted, q_list[i], chosen, noChosenNb, label_q_list[i]);
	
		chosen[sorted[0]] = false;
		added = std::map<int, int>();
		auto skipNodes_t = std::vector<std::set<int> >(sorted.size(), std::set<int>());
		skipNodes = std::vector<std::vector<std::set<int> > >(sorted.size(), skipNodes_t);
		subgraphSearch(nSorted, embedding, tmpcand, label_q_list[i], q_list[i], 0);


		if (embeddingList.size() >= k)
		{
			timer.stop();
			double time = timer.GetRuntime();
			showEmbedding(time, showflag);
			continue;
		}
		else if (embeddingList.size() == 0)
		{
			std::cout << "no matchings0" << std::endl;
			showCount();
			continue;
		}

		//if could not find k non-overlap Embedding, we begin level based search
		candSetList = std::vector<std::set<int> >(sorted.size(), std::set<int>());
		int startIndex = 0;
		while (++overlapSize < sorted.size())
		{
			std::cout << "OverlapSize: " << overlapSize << std::endl;
			std::vector<int> overlapQ(sorted.size(), 0);
			std::vector<std::vector<int> > QoverlapList;
			genOverlapQ(0, overlapSize, overlapQ, QoverlapList);

			genCandSetList3(candSetList, startIndex, label_q_list[i]);

			startIndex = embeddingList.size();


			for (int j = 0; j < QoverlapList.size(); ++j)
			{
				//newSorted re-sort query nodes
				std::vector<std::pair<int, int> > newSorted;
				std::vector<int> &overlapQ = QoverlapList[j];

				chosen = std::vector<bool>(overlapQ.size(), false);

				//store the iterator of data sets for chosen overlap nodes
				for (int m = 0; m < overlapQ.size(); ++m)
				{
					if (overlapQ[m] == 1)
					{
						chosen[sorted[m]] = true;
					}
					else
					{
						int candsize = candidates[label_q_list[i][sorted[m]]].size();
					}
				}

				std::vector<bool> noChosenNb(sorted.size(), false);
				lbnum = labelnumber;
				lbRm = std::vector<int>(sorted.size(), 0);
				nbRm = std::vector<int>(sorted.size(), 0);
				getPairSorted3(rsorted, sorted, newSorted, q_list[i], chosen, noChosenNb, label_q_list[i]);

				std::vector<int> overlapEmb(sorted.size(), -1);
				std::vector<std::vector<int> > overlapEmbList;
				
				if(debugflag)
				{
					std::cout << "newSorted:";
					for (int m = 0; m < newSorted.size(); ++m)
					{
						std::cout << newSorted[m].first << "(" << chosen[newSorted[m].first] << ")(f:" << newSorted[m].second << ")";
					}
					std::cout << std::endl;
				}
		
				getAllEmbList2(0, overlapEmb, overlapEmbList, chosen, newSorted, candSetList, q_list[i], label_q_list[i]);
				
				for (int m = 0; m < overlapEmbList.size(); ++m)
				{
					std::vector<int> & ovlpEmb = overlapEmbList[m];
					added.clear();

					for(auto sn = skipNodes[0].begin(); sn != skipNodes[0].end(); ++sn)
						sn->clear();

					subgraphSearch(newSorted, ovlpEmb, tmpcand, label_q_list[i], q_list[i], overlapSize);
					if (embeddingList.size() >= k)
					{
						timer.stop();
						double time = timer.GetRuntime();
						showEmbedding(time, showflag);
						break;
					}
				}
				if (embeddingList.size() >= k)
					break;
			}
			if (embeddingList.size() >= k)
				break;
		}
		if (embeddingList.size() >= k)
			continue;
		if (!embeddingList.size())
		{
			std::cout << "no matchings" << std::endl;
			showCount();
		}
		else
		{
			timer.stop();
			double time = timer.GetRuntime();
			showEmbedding(time, showflag);
		}
	}
}


//just choose the first overlap node, other nodes are generated in subgraphsearch process
void LevelSub::getAllEmbList2(int startIndex, std::vector<int> &overlapEmb, std::vector<std::vector<int> > &overlapEmbList,
	std::vector<bool> &chosen, std::vector<std::pair<int, int> > &newSorted, std::vector<std::set<int> > &candSetList, Graph &q,
	std::vector<int>& label_q)
{
	int i = startIndex, u;
	for (; i < newSorted.size(); ++i)
	{
		u = newSorted[i].first;
		if (chosen[u] && overlapEmb[u] == -1)
			break;
	}
	if (i >= newSorted.size())
	{
		overlapEmbList.push_back(overlapEmb);
		return;
	}
	for (auto v_t = candSetList[u].begin(); v_t != candSetList[u].end(); ++v_t)
	{
		if (isJoinble(q, u, *v_t, overlapEmb, label_q))
		{
			overlapEmb[u] = *v_t;
			overlapEmbList.push_back(overlapEmb);
			overlapEmb[u] = -1;
		}
	}
}

bool LevelSub::getAllEmbList(int startIndex, std::vector<int> &overlapEmb, std::vector<std::vector<int> > &overlapEmbList, 
	std::vector<bool> &chosen, std::vector<std::pair<int, int> > &newSorted, std::vector<std::set<int> > &candSetList, Graph &q, 
	std::vector<int>& label_q, std::vector<bool> &noChosenNb, int mincand, std::set<int> &added)
{
	int i = startIndex, u, u_f;
	for (; i < newSorted.size(); ++i)
	{
		u = newSorted[i].first;
		if (chosen[u] && overlapEmb[u] == -1)
			break;
	}
	if (i >= newSorted.size())
	{
		overlapEmbList.push_back(overlapEmb);
		return true;
	}

	u_f = newSorted[i].second;
	if (u_f == -1 || overlapEmb[u_f] == -1)
	{
		for (auto v_t = candSetList[label_q[u]].begin(); v_t != candSetList[label_q[u]].end(); ++v_t)
		{
			if (added.find(*v_t) == added.end() && isJoinble(q, u, *v_t, overlapEmb, label_q))
			{
				overlapEmb[u] = *v_t;
				added.insert(*v_t);
				bool flag = getAllEmbList(i + 1, overlapEmb, overlapEmbList, chosen, newSorted, candSetList, q, label_q, noChosenNb, mincand, added);
				overlapEmb[u] = -1;
				added.erase(*v_t);
				//not connected any more, so do not need more
				if (flag && noChosenNb[u])
					return true;
			}
		}
	}
	else
	{
		std::vector<int> &cand = neighborhoodSignature[overlapEmb[u_f]][label_q[u]];
		for (auto v_t = cand.begin(); v_t != cand.end(); ++v_t)
		{
			if (!matched[*v_t] && candSetList[label_q[u]].find(*v_t) == candSetList[label_q[u]].end())
				continue;
			if (added.find(*v_t) == added.end() && isJoinble(q, u, *v_t, overlapEmb, label_q))
			{
				overlapEmb[u] = *v_t;
				added.insert(*v_t);
				bool flag = getAllEmbList(i + 1, overlapEmb, overlapEmbList, chosen, newSorted, candSetList, q, label_q, noChosenNb, mincand, added);
				overlapEmb[u] = -1;
				added.erase(*v_t);
				//not connected any more, so do not need more
				if (flag && noChosenNb[u])
					return true;
			}
		}
	}
	return false;
}

void LevelSub::getPairSorted(std::vector<int> &overlapQ, std::vector<int> &sorted, std::vector<std::pair<int, int> >& newSorted, Graph &q)
{
	std::vector<int> added(overlapQ.size(), -1);
	for (int m = 0; m < overlapQ.size(); ++m)
	{
		if (overlapQ[m] == 1)
		{
			newSorted.push_back(std::pair<int, int>(sorted[m], -1));
			added[sorted[m]] = 1;
			break;
		}
	}
	int m = 0;
	while (newSorted.size() < overlapQ.size())
	{
		int sortedSize = newSorted.size();
		for (; m < sortedSize; ++m)
		{
			int u_f = newSorted[m].first;
			adjacency_iterator_t adj_t, adj_end;
			for (tie(adj_t, adj_end) = adjacent_vertices(u_f, q); adj_t != adj_end; ++adj_t)
			{
				if (added[*adj_t] < 0)
				{
					newSorted.push_back(std::pair<int, int>(*adj_t, u_f));
					added[*adj_t] = 1;
				}
			}
		}
	}
}

void LevelSub::getPairSorted3(std::vector<int> &rsorted, std::vector<int> &sorted, std::vector<std::pair<int, int> >& newSorted,
	Graph &q, std::vector<bool> &chosen, std::vector<bool> &noChosenNb, std::vector<int> &label_q)
{
	int qsize = chosen.size();
	std::vector<int> padded(qsize, -1);

	//sort chosen nodes based on the chosen neighbors and candsets
	for (int m = 0; m < qsize; ++m)
	{
		int u = sorted[m];
		if (chosen[u])
		{
			int countnchsn = 0;
			int countchsn = 0;
			for (auto u_n_t = adjacent_vertices(u, q).first; u_n_t != adjacent_vertices(u, q).second; ++u_n_t)
			{
				if (!chosen[*u_n_t])
					++countnchsn;
				else
					++countchsn;
			}
			if (countnchsn == 0)
				noChosenNb[u] = true;
			else if (newSorted.size() == 0)
			{
				newSorted.push_back(std::pair<int, int>(u, -1));
				lbRm[u] = --lbnum[label_q[u]];
				nbRm[u] = degree(u,q);
				padded[u] = 2;
			}
		}
	}

	if (newSorted.size() == 0)
	{
		for (int i = 0; i < qsize; ++i)
		{
			if (chosen[sorted[i]])
			{
				int u = sorted[i];
				newSorted.push_back(std::pair<int, int>(u, -1));
				lbRm[u] = --lbnum[label_q[u]];
				nbRm[u] = degree(u,q);
				padded[u] = 2;
				break;
			}
		}
	}

	std::vector<std::pair<int, int> > noChosenNbList;

	int m = 0;
	while (newSorted.size() <qsize)
	{
		int sortedSize = newSorted.size();
		for (; m < sortedSize; ++m)
		{
			int u_f = newSorted[m].first;
			std::vector<int> nb;
			adjacency_iterator_t adj_t, adj_end;
			for (tie(adj_t, adj_end) = adjacent_vertices(u_f, q); adj_t != adj_end; ++adj_t)
			{
				if (padded[*adj_t] < 0)
				{
					if (noChosenNb[*adj_t])
						noChosenNbList.push_back(std::pair<int, int>(*adj_t, u_f));
					else
						nb.push_back(*adj_t);
					padded[*adj_t] = 1;
				}
			}
			std::sort(nb.begin(), nb.end(), [&rsorted](int a, int b){return rsorted[a] < rsorted[b]; });
			for (int i = 0; i < nb.size(); ++i)
			{
				int u = nb[i];
				newSorted.push_back(std::pair<int, int>(u, u_f));
				lbRm[u] = --lbnum[label_q[u]];
				for (tie(adj_t, adj_end) = adjacent_vertices(u, q); adj_t != adj_end; ++adj_t)
				{
					if(padded[*adj_t] == 1 || padded[*adj_t] < 0)
						++nbRm[u];
				}
				padded[u] = 2;
			}
		}
		if (m == newSorted.size())
		{
			for (int i = noChosenNbList.size() - 1; i >= 0; --i)
			{
				newSorted.push_back(noChosenNbList[i]);
				int u = noChosenNbList[i].first;
				adjacency_iterator_t adj_t, adj_end;
				for (tie(adj_t, adj_end) = adjacent_vertices(u, q); adj_t != adj_end; ++adj_t)
				{
					if(padded[*adj_t] == 1 || padded[*adj_t] < 0)
						++nbRm[u];
				}
				lbRm[u] = --lbnum[label_q[u]];
				padded[u] = 2;
			}
			noChosenNbList.clear();
		}
	}
}

void LevelSub::getPairSorted2(std::vector<int> &overlapQ, std::vector<int> &sorted, std::vector<std::pair<int, int> >& newSorted, Graph &q)
{
	std::vector<int> added(overlapQ.size(), -1);
	bool first = true;
	for (int m = 0; m < overlapQ.size(); ++m)
	{
		if (overlapQ[m] == 1)
		{
			added[sorted[m]] = 1;
			if (first)//just put the first overlap node into newsorted
			{
				newSorted.push_back(std::pair<int, int>(sorted[m], -1));
				added[sorted[m]] = 0;
				first = false;
			}
		}
	}
	int m = 0;
	while (newSorted.size() < overlapQ.size())
	{
		int sortedSize = newSorted.size();
		for (; m < sortedSize; ++m)
		{
			int u_f = newSorted[m].first;
			adjacency_iterator_t adj_t, adj_end;
			for (tie(adj_t, adj_end) = adjacent_vertices(u_f, q); adj_t != adj_end; ++adj_t)
			{
				if (added[*adj_t] < 0)
				{
					newSorted.push_back(std::pair<int, int>(*adj_t, u_f));
					added[*adj_t] = 0;
				}
				else if (added[*adj_t] == 1)//overlap node
				{
					newSorted.push_back(std::pair<int, int>(*adj_t, -1));
					added[*adj_t] = 0;
				}
			}
		}
	}
}


void LevelSub::genOverlapQ(int start, int remain, std::vector<int> &overlapQ, std::vector<std::vector<int> > &QoverlapList)
{
	if (remain == 0)
	{
		QoverlapList.push_back(overlapQ);
		return;
	}
	if (remain + start > overlapQ.size())
		return;
	overlapQ[start] = 1;
	genOverlapQ(start + 1, remain - 1, overlapQ, QoverlapList);
	overlapQ[start] = 0;
	genOverlapQ(start + 1, remain, overlapQ, QoverlapList);
}

void LevelSub::genCandSetList(std::vector<std::set<int> > &candSetList, int startIndex, std::vector<int> &label_q)
{
	for (int i = startIndex; i < embeddingList.size(); ++i)
	{
		for (int j = 0; j < embeddingList[i].size(); ++j)
		{
			//candSetList[j].insert(embeddingList[i][j]);
			candSetList[label_q[j]].insert(embeddingList[i][j]);
		}
	}
}

void LevelSub::genCandSetList3(std::vector<std::set<int> > &candSetList, int startIndex, std::vector<int> &label_q)
{
	if(debugflag)
		std::cout << "test in" << std::endl;
	for (int i = startIndex; i < embeddingList.size(); ++i)
	{
		std::map<int, std::set<int> > lCandidate;
		for (int j = 0; j < embeddingList[i].size(); ++j)
		{
			//candSetList[j].insert(embeddingList[i][j]);
			lCandidate[label_q[j]].insert(embeddingList[i][j]);
		}
		for (int j = 0; j < embeddingList[i].size(); ++j)
		{
			//candSetList[j].insert(embeddingList[i][j]);
			for(auto l_t = lCandidate[label_q[j]].begin(); l_t != lCandidate[label_q[j]].end(); ++l_t)
			{
				if(nsFilter(*l_t, j))
					candSetList[j].insert(*l_t);
			}
		}
	}
	if(debugflag)
		std::cout << "test out" << std::endl;
}

void LevelSub::genCandSetList2(std::vector<std::set<int> > &candSetList, int startIndex)
{
	for (int i = startIndex; i < embeddingList.size(); ++i)
	{
		for (int j = 0; j < embeddingList[i].size(); ++j)
		{
			//candSetList[j].insert(embeddingList[i][j]);
			candSetList[j].insert(embeddingList[i][j]);
		}
	}
}


void LevelSub::subgraphSearch(std::vector<std::pair<int, int> > &sorted, std::vector<int> &embedding, std::vector<std::vector<int> > &cand,
	std::vector<int> &label_q, Graph &q, int overlapSize)
{
	++nodeGraphSearchCount;


	int i = 0; 
	while (i < sorted.size() && embedding[sorted[i].first] != -1)
	{
		added[embedding[sorted[i].first]] = sorted[i].first;
		++i;
	}
	if (i == sorted.size())
	{
		return;
	}

	for(auto t = skipNodes[i].begin(); t != skipNodes[i].end(); ++t)
		t->clear();

	int u = sorted[i].first;
	int u_f = sorted[i].second;
	int label = label_q[u];

	std::vector<bool> conflict(gCT[u]);

	std::vector<int>::iterator v_t, v_end;
	bool candFlag = false;
	if (u_f != -1 && embedding[u_f] != -1)
	{
		v_t = neighborhoodSignature[embedding[u_f]][label].begin();
		v_end = neighborhoodSignature[embedding[u_f]][label].end();
	}
	else
	{
		candFlag = true;
		v_t = cand[label].begin();
		v_end = cand[label].end();
	}

	if(debugflag)
	{
		std::cout << "subgraphSearch: Searching matches for u" << u;
	}

	bool nextlayer = false;
	bool hasSkip = false;
	int validCandidateCount = 0;
	for (; v_t != v_end; ++v_t)
	{

		if(*v_t == -1)
			continue;
		else if(chosen[u] && (!matched[*v_t] || candSetList[u].find(*v_t) == candSetList[u].end()) )
			continue;
		else if(!chosen[u] && matched[*v_t])
		{
			if (added.find(*v_t) != added.end())
			{
				conflict[added[*v_t]] = true;
			}
			continue;
		}

		if(debugflag)
			std::cout << "matching " << *v_t << " for " << u << std::endl;

		if (added.find(*v_t) != added.end())
		{
			conflict[added[*v_t]] = true;
			continue;
		}


		if (!chosen[u])
		{
			if (degree(*v_t, g) < degree(u, q) || !nsFilter(*v_t, u))
			{
				if (candFlag)
					*v_t = -1;
				continue;
			}
		}

		
		bool skipFlag = false;
		for(int j = 0; j < i; ++j)
		{
			if(skipNodes[j][u].find(*v_t) != skipNodes[j][u].end())
			{
				if(debugflag)
					std::cout << "skipnode" << std::endl;
				skipFlag = true;
				break;
			}
		}
		if(skipFlag)
		{
			hasSkip = true;
			continue;
		}


		if (isJoinble(q, u, *v_t, embedding, label_q))
		{
			nextlayer = true;

			if(nbRm[u] == 0 && lbRm[u] == validCandidateCount && !candFlag && chosen[u])
				v_end = v_t + 1;
			++validCandidateCount;

			std::vector<int> oldEmbedding = embedding;
			std::map<int, int> oldadded = added;
			embedding[u] = *v_t;
			matched[*v_t] = 1;
			added[*v_t] = u;

			if (chosen[u])
			{
				subgraphSearch(sorted, embedding, cand, label_q, q, overlapSize);
				if (embeddingList.size() >= k)
					return;
			}
			else if (subgraphSearchD(i+1, sorted, embedding, cand, label_q, q))
			{
				falseu = -1;
				gconflict.clear();

				embeddingList.push_back(embedding);
				showEmbedding(embedding);

				if (embeddingList.size() >= k)
					return;

				embedding = oldEmbedding;
				added = oldadded;
				if (!chosen[u])
					continue;
			}

			embedding[u] = -1;
			added.erase(*v_t);

			if (!chosen[u])
				matched[*v_t] = 0;

			if(nbRm[u] == 0 && lbRm[u] == 0)
				return;

			if(falseu >= i+1)
			{
				int fu = sorted[falseu].first;
				if(!gconflict[u])
				{
					return;
				}
				else if(sorted[falseu].second == u_f && degree(fu, q) == 1 && degree(u, q) == 1)
				{
					if(debugflag)
					{
						std::cout << "falseu:" << sorted[falseu].first << ", u_f:" << sorted[falseu].second << std::endl;
					}
					return;
				}

				if(i != 0 && !gconflict[sorted[i-1].first])
					skipNodes[i-1][u].insert(*v_t);
			}
			if(badv && i != 0)
			{
				skipNodes[i-1][u].insert(*v_t);
			}
			badv = false;
			falseu = -1;
			gconflict.clear();
		}
	}
	if(hasSkip && !nextlayer)
	{
		badv = true;
		nextlayer = true;
	}
	if(!candFlag && !nextlayer)
	{
		falseu = i;
		gconflict = conflict;
		badv = false;
	}
}

void LevelSub::update2(std::vector<int> &embedding, std::vector<int> &label_q)
{
	embeddingList.push_back(embedding);
}

bool LevelSub::subgraphSearchD(int i, std::vector<std::pair<int, int> > &sorted, std::vector<int> &embedding,
	std::vector<std::vector<int> > &cand, std::vector<int> &label_q, Graph &q)
{
	++nodeGraphSearchCount;
	if (i >= sorted.size())
	{
		return true;
	}

	for(auto t = skipNodes[i].begin(); t != skipNodes[i].end(); ++t)
		t->clear();

	int u = sorted[i].first;
	int u_f = sorted[i].second;
	int label = label_q[u];


	if (embedding[u] != -1)
	{
		added[embedding[u]] = u;
		return subgraphSearchD(i + 1, sorted, embedding, cand, label_q, q);
	}

	std::vector<bool> conflict(gCT[u]);
	
	bool candFlag = false;
	std::vector<int>::iterator v_t, v_end;
	if (u_f != -1 && embedding[u_f] != -1)
	{
		v_t = neighborhoodSignature[embedding[u_f]][label].begin();
		v_end = neighborhoodSignature[embedding[u_f]][label].end();
	}
	else
	{
		candFlag = true;
		v_t = cand[label].begin();
		v_end = cand[label].end();
	}

	if(debugflag)
	{
		std::cout << "subgraphSearchD: Searching matches for u" << u;
		showEmbedding(embedding);
	}
	
	bool nextlevel = false;
	int validCandidateCount = 0;

	for (; v_t != v_end; ++v_t)
	{
		if(*v_t == -1)
			continue;
		else if(chosen[u] && (!matched[*v_t] || candSetList[u].find(*v_t) == candSetList[u].end()) )
			continue;
		else if(!chosen[u] && matched[*v_t])
		{
			if (added.find(*v_t) != added.end())
			{
				conflict[added[*v_t]] = true;
			}
			continue;
		}

		if(debugflag)
			std::cout << "D: matching " << *v_t << " for " << u << std::endl;

		if (added.find(*v_t) != added.end())
		{
			conflict[added[*v_t]] = true;
			continue;
		}

		if (!chosen[u])
		{
			if (degree(*v_t, g) < degree(u, q) || !nsFilter(*v_t, u))
			{
				if (candFlag)
					*v_t = -1;
				continue;
			}
		}
		bool skipFlag = false;
		for(int j = 0; j < i; ++j)
		{
			if(skipNodes[j][u].find(*v_t) != skipNodes[j][u].end())
			{
				if(debugflag)
					std::cout << "skipnode" << std::endl;
				skipFlag = true;
				break;
			}
		}
		if(skipFlag)
		{
			nextlevel = true;
			continue;
		}


		if (isJoinble(q, u, *v_t, embedding, label_q))
		{
			nextlevel = true;
			if(nbRm[u] == 0 && lbRm[u] == validCandidateCount)
				v_end = v_t + 1;
			++validCandidateCount;

			embedding[u] = *v_t;
			added[*v_t] = u;
			matched[*v_t] = 1;
			if (subgraphSearchD(i + 1, sorted, embedding, cand, label_q, q))
				return true;
			embedding[u] = -1;
			added.erase(*v_t);
			if (!chosen[u])
				matched[*v_t] = 0;
			
			
			if(nbRm[u] == 0 && lbRm[u] == 0)
				return false;
			
			//return to it's father or same label one
			if(debugflag)
				std::cout << "falseu" << falseu << " i:" << i << std::endl;
			if(falseu >= i+1)
			{
				int fu = sorted[falseu].first;

				if(!gconflict[u])
				{
					if(debugflag)
					{
						std::cout << "isEdge:" << fu << " " <<  u << ":"<< isEdgeq(fu,label_q[fu], u, label, q) << std::endl;
						std::cout << "label fu:" << label_q[fu] << " label u" <<  label << std::endl;
					}
					return false;
				}
				else if(sorted[falseu].second == u_f && degree(fu, q) == 1 && degree(u, q) == 1)
				{
					if(debugflag)
					{
						std::cout << "falseu:" << sorted[falseu].first << ", u_f:" << sorted[falseu].second <<  std::endl;
					}
					return false;
				}

				if(i != 0 && !gconflict[sorted[i-1].first])
					skipNodes[i-1][u].insert(*v_t);
			}
			falseu = -1;
			gconflict.clear();
		}
	}
	if(debugflag)
	{
		std::cout << "candFlag:" << candFlag << "Nextlevel:" << nextlevel << std::endl;
	}

	if(!candFlag && !nextlevel)
	{
		falseu = i;
		gconflict = conflict;
	}
	return false;
}

inline void LevelSub::showEmbedding(double time, bool flag)
{
	std::cout << "There are " << embeddingList.size() << " embeddings" << std::endl;
	std::vector<std::set<int> > candSetList(embeddingList[0].size(), std::set<int>());
	genCandSetList2(candSetList, 0);
	std::set<int> coverage;

	for (int i = 0; i < candSetList.size(); ++i)
	{
		std::cout << i << " " << candSetList[i].size() << std::endl;
		std::cout << "[";
		for (auto s_t = candSetList[i].begin(); s_t != candSetList[i].end(); ++s_t)
		{
			std::cout << *s_t << " ";
			coverage.insert(*s_t);
		}
		std::cout << "]" << std::endl;
	}


	if (flag)
	{
		for (int i = 0; i < embeddingList.size(); ++i)
		{
			std::cout << "{";
			for (int j = 0; j < embeddingList[i].size(); ++j)
				std::cout << j << "<>" << embeddingList[i][j] << " ";
			std::cout << "}" << std::endl;
		}
	}
	showCount();
	std::cout << "Embedding Size:" << embeddingList.size() << ", q size:" << chosen.size() << std::endl;
	std::cout << "Cover nodes:" << coverage.size() << ", Coverage(#nodes/q_size*embeddings):" << ((double)(coverage.size()))/(chosen.size()*embeddingList.size()) << std::endl;
	std::cout << "Cost time(ms) " << time  << std::endl;
}

inline void LevelSub::showEmbedding(const std::vector<int> &embedding)
{
	std::cout << "{";
	for (int j = 0; j < embedding.size(); ++j)
		std::cout << j << "<>" << embedding[j] << " ";
	std::cout << "}" << std::endl;
}

//neighborhood signature
bool LevelSub::nsFilter(vertex_t v, vertex_t u)
{
	++nsFilterCount;
	std::map<int, std::set<int> >::iterator ns_qt;
	for (ns_qt = nsQ[u].begin(); ns_qt != nsQ[u].end(); ++ns_qt)
	{
		auto ns_gt = neighborhoodSignature[v].find(ns_qt->first);
		if (ns_gt == neighborhoodSignature[v].end() || ns_gt->second.size() < ns_qt->second.size())
		{
			return false;
		}
	}
	return true;
}

//test if the vertex can be join to the result
bool LevelSub::isJoinble(Graph &q, vertex_t u, vertex_t v, std::vector<int> &embedding, std::vector<int> &label_q)
{
	++isJoinbleCount;
	adjacency_iterator_t adj, adj_end;

	for (tie(adj, adj_end) = adjacent_vertices(u,q); adj != adj_end; ++adj)
	{
		if (embedding[*adj] != -1)
		{
			if (!isEdge(v, label_q[u], embedding[*adj], label_q[*adj]))
				return false;
		}
	}
	return true;
}

//maybe need to find a more efficient way to test edge
//need to study the graph library more to choose the best container
bool LevelSub::isEdge(vertex_t v, int lv, vertex_t u, int lu)
{
	isEdgeCount++;

	int label = lu;


//	if (degree(v, g) > degree(u, g))
	if (neighborhoodSignature[v][lu].size() > neighborhoodSignature[u][lv].size())
	{
		vertex_t tmp = v;
		v = u;
		u = tmp;
		label = lv;
	}

	for(auto adj = neighborhoodSignature[v][label].begin(); adj != neighborhoodSignature[v][label].end(); ++adj)
	{
		if (u == *adj)
			return true;
	}
	return false;
}

bool LevelSub::isEdgeq(vertex_t v, int lv, vertex_t u, int lu, Graph &q)
{
	isEdgeCount++;


	adjacency_iterator_t adj, adj_end;

	if (degree(v, q) > degree(u, q))
	{
		vertex_t tmp = v;
		v = u;
		u = tmp;
	}

	for (tie(adj, adj_end) = adjacent_vertices(v, q); adj != adj_end; ++adj)
	{
		if (u == *adj)
			return true;
	}
	return false;
}

void LevelSub::showCount()
{
	std::cout << "isJoinbleCount:" << isJoinbleCount << std::endl;
	std::cout << "GraphSearchCount:" << nodeGraphSearchCount << std::endl;
	std::cout << "nsFilterCount:" << nsFilterCount << std::endl;
	std::cout << "isEdgeCount:" << isEdgeCount << std::endl;
	isJoinbleCount = nodeGraphSearchCount = nsFilterCount = isEdgeCount = 0;
}
