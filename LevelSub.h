#ifndef LEVEL_SUB
#define LEVEL_SUB
#include <string>
#include <vector>
#include <map>
#include <queue>

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/next_prior.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/dynamic_bitset.hpp>

using namespace boost;

/*
struct label_t{
typedef vertex_property_tag kind;
};


struct Vertex{
int id;
int label;
};
*/
//typedef property<label_t, int> LabelProperty;
typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int>, property<edge_index_t, int> > Graph;
//typedef adjacency_list<vecS, vecS, undirectedS, Vertex, property<edge_index_t, int> > Graph;
typedef graph_traits<Graph>::vertex_descriptor vertex_t;
typedef graph_traits<Graph>::edge_descriptor edge_t;
typedef graph_traits<Graph>::vertex_iterator vertex_iterator_t;
typedef graph_traits<Graph>::edge_iterator edge_iterator_t;
typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator_t;

typedef dynamic_bitset<> d_bitset;

/* KJH */
typedef adjacency_list<vecS, vecS, directedS, property<vertex_index_t, int>, property<edge_index_t, int> > dirGraph;
typedef graph_traits<dirGraph>::vertex_descriptor dir_vertex_t;
typedef graph_traits<dirGraph>::edge_descriptor dir_edge_t;
typedef graph_traits<dirGraph>::vertex_iterator dir_vertex_iterator_t;
typedef graph_traits<dirGraph>::edge_iterator dir_edge_iterator_t;
typedef graph_traits<dirGraph>::adjacency_iterator dir_adjacency_iterator_t;

typedef std::map<dir_vertex_t, std::vector<vertex_t> > NEC_vertex_map_t;
typedef std::map<std::pair<dir_vertex_t,vertex_t>, std::vector<vertex_t> > Cand_sub_t;

typedef std::map<std::pair<dir_vertex_t,vertex_t>, std::vector<vertex_t> >::iterator my_it;

class LevelSub{
	typedef std::set<int>::iterator set_t;
private:
	std::string DataGraphFileName;
	std::string QueryGraphFileName;
	int k;

	//std::vector<int>* label_q_p;
	std::vector<Graph> q_list;//queries 

	/* KJH */
	std::vector<dirGraph> NEC_tree_list; 
	std::vector<NEC_vertex_map_t> NEC_data_list;
	std::vector<dir_vertex_t> NEC_vertices;
	std::map<dir_vertex_t,dir_vertex_t> parent_map;
	Cand_sub_t cand_sub;
	std::map<dir_vertex_t,std::set<vertex_t> > cand_key_pair;
	std::map<vertex_t,bool> data_visited;
	std::map<dir_vertex_t,int> NEC_label;
	std::map<vertex_t,int> label_v;
	std::vector<std::vector<dir_vertex_t> > matching_order_list;
	std::map<dir_vertex_t, std::vector<int> > curr_comb_idx_map;
	std::map<dir_vertex_t, std::vector<vertex_t> > curr_comb_map;
	std::map<dir_vertex_t, std::vector<int> > curr_perm_idx_map;
	std::map<vertex_t, vertex_t> map_table;

	std::vector<std::vector<int> > label_q_list;
	std::vector<std::map<int, std::set<int> > > nsQ;
	std::vector<std::vector<int> > embeddingList;

	std::vector<bool> chosen;
	std::map<int, int> lbnum;
	std::vector<int> lbRm; // label remain
	std::vector<int> nbRm;// neighbor remain
	std::vector<std::set<int> > candSetList;
	std::map<int, int> added;

	int falseu = -1;
	std::vector<bool>  gconflict;
	std::vector<std::vector<bool> > gCT; //global conflict table
	//int gMaxConflict = 0;

	//std::vector<std::set<int> > skipNodes;
	std::vector<std::vector<std::set<int> > > skipNodes;
	int badv = false;


	bool showflag = false;
	bool debugflag = false;

	Graph g;//Data graph
	std::vector<std::map<int, std::vector<int> > > neighborhoodSignature;
	std::vector<std::vector<int> > candidates;
	d_bitset matched;//global matched mask => no overlap

	long isJoinbleCount = 0;
	long nodeGraphSearchCount = 0;
	long nsFilterCount = 0;
	long isEdgeCount = 0;

	void LoadGraph();
	void LoadQuery();
	void buildSignature(int v_label, int v, int u);
	void randomShuffle();
	void queryOrder(Graph &q, std::vector<int> &label_q, std::vector<int>& sorted);
	void genOverlapQ(int start, int remain, std::vector<int> &overlapQ, std::vector<std::vector<int> > &QoverlapList);

	//void subgraphSearch(std::vector<int> &sorted, std::vector<int> &embedding, std::vector<std::vector<int> > &tmpcand, 
		//std::vector<int> &label_q, Graph &q);
	void subgraphSearch(std::vector<std::pair<int, int> > &sorted, std::vector<int> &embedding, std::vector<std::vector<int> > &tmpcand,
		std::vector<int> &label_q, Graph &q, int overlapSize);

	bool getAllEmbList(int startIndex, std::vector<int> &overlapEmb, std::vector<std::vector<int> > &overlapEmbList, 
		std::vector<bool> &chosen,std::vector<std::pair<int, int> > &newSorted, std::vector<std::set<int> > &candSetList,
		Graph &q, std::vector<int> &label_q, std::vector<bool> &noChosenNb, int mincand, std::set<int> &added);
	void getAllEmbList2(int startIndex, std::vector<int> &overlapEmb, std::vector<std::vector<int> > &overlapEmbList,
		std::vector<bool> &chosen, std::vector<std::pair<int, int> > &newSorted, std::vector<std::set<int> > &candSetList,
		Graph &q, std::vector<int> &label_q);

	void getPairSorted(std::vector<int> &overlapQ, std::vector<int> &sorted, std::vector<std::pair<int, int> >& newSorted, Graph &q);
	void getPairSorted2(std::vector<int> &overlapQ, std::vector<int> &sorted, std::vector<std::pair<int, int> >& newSorted, Graph &q);
	void getPairSorted3(std::vector<int> &rsorted, std::vector<int> &sorted, 
		std::vector<std::pair<int, int> >& newSorted, Graph &q, std::vector<bool> &chosen, std::vector<bool> &noChosenNb, std::vector<int> &label_q/*, std::vector<bool> &allChosenNb*/);

	void update2(std::vector<int> &embedding, std::vector<int> &label_q);

	bool subgraphSearchD(int i, std::vector<std::pair<int, int> >& sorted, std::vector<int> &embedding,
		std::vector<std::vector<int> > &cand, std::vector<int> &label_q, Graph &q);
	//bool subgraphSearchD(int i, std::vector<int> &sorted, std::vector<int> &embedding,
		//std::vector<std::vector<int> > &cand, std::vector<int> &label_q, Graph &q);
	void genCandSetList(std::vector<std::set<int> > &candSetList, int startIndex, std::vector<int> &label_q);
	void genCandSetList3(std::vector<std::set<int> > &candSetList, int startIndex, std::vector<int> &label_q);
	void genCandSetList2(std::vector<std::set<int> > &candSetList, int startIndex);

	bool nsFilter(vertex_t v, vertex_t u);
	bool isJoinble(Graph &q, vertex_t u, vertex_t v, std::vector<int> &embedding, std::vector<int> &label_q);
	bool isEdge(vertex_t v, int lv, vertex_t u, int lu);
	bool isEdgeq(vertex_t v, int lv, vertex_t u, int lu, Graph &q);
public:
	LevelSub(const std::string DataFileName, const std::string QueryFileName, const int k, const int debug);
	void subgraphSearch();
	void showEmbedding(double time, bool flag);
	void showEmbedding(const std::vector<int> &embedding);
	void showCount();
	~LevelSub();

	void rewriteToNECTree(); // KJH : Rewrite to NEC Tree
	std::vector<std::vector<vertex_t> > FindNEC(std::vector<vertex_t> v);
	bool NLF(vertex_t u, vertex_t v);
	bool exploreCR(dir_vertex_t u, std::vector<vertex_t> V, vertex_t v_parent);
	void clearCR(dir_vertex_t u, vertex_t v);
	int cand_sub_size(dir_vertex_t u);
	double non_tree_edge_num(dir_vertex_t u);
	std::vector<dir_vertex_t> find_path(dir_vertex_t u);
	void matchingOrder();
	std::vector<vertex_t> NextComb(dir_vertex_t u,std::vector<vertex_t> C);
	void GenPerm(int idx);
	bool NextPerm(dir_vertex_t u_i);
	bool IsJoinable(vertex_t v, vertex_t C_i);
	void SubgraphSerach();
};
#endif
