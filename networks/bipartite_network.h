#ifndef BIPARTITENETWORK_H
#define BIPARTITENETWORK_H

#include "weighted_network.h"
#include <queue>

using namespace std; 

class BipartiteNetwork : public WeightedNetwork {
	public:
		

		int B; //indexing will be left node:(0,1, ..., B-1) right nodes:(B, B+1, ..., num_nodes-1)
		map<int, int> colour;	//colour vector for bipartite check
					
		BipartiteNetwork():WeightedNetwork(){
			network_type="bipartite_network";
			B = 0;
		}		
		
		
		BipartiteNetwork(const char* filename);
		bool check_bipartite(map< int, vector<link> > &lk);
		bool is_bipartite();
		void compute_indices(map< int, vector<link> > &lk);
		BipartiteNetwork(const WeightedNetwork &G );
		BipartiteNetwork induced_graph_side(vector<int> &label, bool left);
		void print_basic(bool use_original_ids=true);

		
};


BipartiteNetwork::BipartiteNetwork(const WeightedNetwork &G) : WeightedNetwork(G){
	network_type="bipartite_network";	
}

BipartiteNetwork::BipartiteNetwork(const char* filename){
	network_type="bipartite_network";	
	map< int, vector<link> > lk = read(filename);
	compute_indices(lk);
	reindex(lk);
}

//print some stuff
void BipartiteNetwork::print_basic(bool use_original_ids){

	Network::print_basic(use_original_ids);
	cout << "Number of left nodes = " << B << endl;
	cout << "Number of right nodes = " << num_nodes - B << endl;
}


///See if the network is bipartite
bool BipartiteNetwork::is_bipartite(){
	map< int, vector<link> > lk = to_map();
	return check_bipartite(lk);
}

///See if a map is bipartite
bool BipartiteNetwork::check_bipartite(map< int, vector<link> > &lk){
	
	//colouring algorithm, will not work for disconnected networks
	queue<int> q; q.push( lk.begin()->first );
	colour.clear(); colour[ lk.begin()->first ] = 1;

	while(!q.empty()){
		int i = q.front(); q.pop();
		for(auto &j : lk[i]){ //int j=0; j<lk[i].size(); ++j){
			int nn = j.first; //the neighbour node
			if(colour.find(nn) != colour.end()){ //have we already coloured it?
				if(colour[i] == colour[nn]){ return false; }
			} else {
				colour[nn] = 1-colour[i];
				q.push(nn);
			}
		}
	}
	return true;

}

///compute indices of bipartite network.
///put all left nodes at [0,1,...,B-1]
///right nodes at [B, B+1, ..., num_nodes-1]
void BipartiteNetwork::compute_indices(map< int, vector<link> > &lk){
	
	bool bip = check_bipartite(lk);

	//put original first id on left!
	int left_colour = colour[ lk.begin()->first ]; 

	B = 0;
	new_to_old.resize(num_nodes);		
	for(map<int, int>::iterator it=colour.begin(); it != colour.end(); ++it){ 
		if( it->second == left_colour ){ new_to_old[B] = it->first;  old_to_new[it->first] = B; ++B; }
	}
	int C = B;
	for(map<int, int>::iterator it=colour.begin(); it != colour.end(); ++it){ 
		if( it->second != left_colour ){ new_to_old[C] = it->first;  old_to_new[it->first] = C; ++C; }
	}		

}

BipartiteNetwork BipartiteNetwork::induced_graph_side(vector<int> &label, bool left){
	
	if(label.size() != (unsigned)num_nodes){
		if(left){
			if(unsigned(B) != label.size()){
				cerr << "Trying to induce " << B << " nodes with " << label.size() << " labels" << endl;
				exit(1);
			}
		} else {
			if(unsigned(num_nodes - B) != label.size()){
				cerr << "Trying to induce " << num_nodes - B << " nodes with " << label.size() << " labels" << endl;
				exit(1);
			}
		}
	} 
	
	int ll = (left) ? 0 : B;
	int ul = (left) ? B : num_nodes; 
	int lab_idx = (left) ? 0 : num_nodes;
	
	map< int, map<int, double> > mp;
	for(int i=ll; i<ul; ++i){
		int a = lab_idx + label[ i-ll ]; 
		for(auto &j : get_nbrs[i]){ 
			int b = j.first; 
			if(mp.find(a) == mp.end()){ mp[a][b] = j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = j.second;  }
				else{ mp[a][b] += j.second; }
			}
		}
	}

	
	map< int, vector<link> > lk; 
	for( auto it = mp.begin(); it != mp.end(); ++it){
		for(auto jt = it->second.begin(); jt != it->second.end(); ++jt){

			int a = it->first;
			int b = jt->first;
			int w = jt->second;
			if(lk.find(a) == lk.end()){
				lk[ a ] = vector<link>{ make_pair(b, w) };
			} else {
				lk[ a ].push_back( make_pair(b, w) );
			}
			
			if(lk.find(b) == lk.end()){
				lk[ b ] = vector<link>{ make_pair(a, w) };
			} else {
				lk[ b ].push_back( make_pair(a, w) );
			}	

		}
	}

	
	BipartiteNetwork bp = BipartiteNetwork(lk);
	bp.compute_indices(lk);
	bp.reindex(lk);
	
	return bp;
}


#endif
