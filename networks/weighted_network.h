#ifndef WEIGHTEDNETWORK_H
#define WEIGHTEDNETWORK_H

#include "network.h"
#include <unordered_map>
#include <unordered_set>

using namespace std; 

///An undirected weighted network
class WeightedNetwork : public Network {
	public:
				
		WeightedNetwork():Network(){
			network_type="weighted_network";
		}		
		
		WeightedNetwork(const char* filename);
		WeightedNetwork(adjacency_map &lk);
		WeightedNetwork(const Network &G );
	
		WeightedNetwork induced_graph(unordered_map<int, int> &label);
};



///Use provided node labels to compute induced graph
WeightedNetwork WeightedNetwork::induced_graph(unordered_map<int, int> &label){ 
	
	if( (unsigned)num_nodes != label.size()){
		cerr << "Trying to induce a " << num_nodes << " network with " << label.size() << " labels" << endl;
		exit(1);
	}
	
	adjacency_map mp;
	for(auto &i : label){
		int a = i.second; 
		for(auto &j : get_nbrs[i.first]){
			int b = label[ j.first ];
			if(mp.find(a) == mp.end()){ mp[a][b] = j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = j.second;  }
				else{ mp[a][b] += j.second; }
			}
		}
	}
	return WeightedNetwork(mp);
	
	
}


WeightedNetwork::WeightedNetwork(const char* filename){
	network_type="weighted_network";	
	read(filename);
}

WeightedNetwork::WeightedNetwork(adjacency_map &mp) : Network(mp){
	network_type="weighted_network";	
}


WeightedNetwork::WeightedNetwork(const Network &G) : Network(G){
	network_type="weighted_network";	
}




#endif
