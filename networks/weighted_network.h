#ifndef WEIGHTEDNETWORK_H
#define WEIGHTEDNETWORK_H

#include "network.h"

using namespace std; 

///An undirected weighted network
class WeightedNetwork : public Network {
	public:
				
		WeightedNetwork():Network(){
			network_type="weighted_network";
		}		
		
		WeightedNetwork(const char* filename);
		WeightedNetwork(map< int, vector<link> > &lk);
		WeightedNetwork(const Network &G );
	
		WeightedNetwork induced_graph(vector<int> &label);
};

///Use provided node labels to compute induces graph
WeightedNetwork WeightedNetwork::induced_graph(vector<int> &label){ 
	if( (unsigned)num_nodes != label.size()){
		cerr << "Trying to induce a " << num_nodes << " network with " << label.size() << " labels" << endl;
		exit(1);
	}
	//Use map because labels can be anything : 1, 1, 1, 3, 4, 5
	//map< map > makes it easier to sum up without having to search a vector each time.
	map< int, map<int, double> > mp; 
	for(int i=0; i<num_nodes; ++i){
		int a = label[ i ]; //label of this node (using the new node ids!)
		for(auto &j : get_nbrs[i]){
			int b = label[ j.first ];
			if(mp.find(a) == mp.end()){ mp[a][b] = j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = j.second;  }
				else{ mp[a][b] += j.second; }
			}
		}
	}

	//use this to set up a new Network
	map< int, vector<link> > lk; 
	for( auto it = mp.begin(); it != mp.end(); ++it){
		auto jt = it->second.begin();
		lk[ it->first ] = vector<link>{ make_pair(jt->first, jt->second) };
		++jt;
		for( ;jt != it->second.end(); ++jt){
			lk[ it->first ].push_back( make_pair(jt->first, jt->second) );
		} 
	}
	
	return WeightedNetwork(lk);
	
}


WeightedNetwork::WeightedNetwork(const char* filename){
	network_type="weighted_network";	
	map< int, vector<link> > lk = read(filename);
	compute_indices(lk);
	reindex(lk);
	
}

WeightedNetwork::WeightedNetwork(map< int, vector<link> > &lk) : Network(lk){
	network_type="weighted_network";	
}


WeightedNetwork::WeightedNetwork(const Network &G) : Network(G){
	network_type="weighted_network";	
}




#endif
