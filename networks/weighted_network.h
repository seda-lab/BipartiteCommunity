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
		int a = i.second; //label[ i ]; //label of this node (using the new node ids!)
		//cout << "induced_graph:: label " << a << " for node " << i << endl;
		for(auto &j : get_nbrs[i.first]){
			int b = label[ j.first ];
			//cout << "induced_graph:: label " << b << " for node " << j.first << endl;
			if(mp.find(a) == mp.end()){ mp[a][b] = j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = j.second;  }
				else{ mp[a][b] += j.second; }
			}
		}
	}
	return WeightedNetwork(mp);
	
	/*//Use map because labels can be anything : 1, 1, 1, 3, 4, 5
	//map< map > makes it easier to sum up without having to search a vector each time.
	map< int, map<int, double> > mp; 
	for(int i=0; i<num_nodes; ++i){
		int a = label[ i ]; //label of this node (using the new node ids!)
		//cout << "induced_graph:: label " << a << " for node " << i << endl;
		for(auto &j : get_nbrs[i]){
			int b = label[ j.first ];
			//cout << "induced_graph:: label " << b << " for node " << j.first << endl;

			if(mp.find(a) == mp.end()){ mp[a][b] = j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = j.second;  }
				else{ mp[a][b] += j.second; }
			}
		}
	}*/
		/*for( auto it = mp.begin(); it != mp.end(); ++it){
		cout << "induced_graph:: community[" << it->first << "] connected to ( ";
		for( auto jt = it->second.begin();jt != it->second.end(); ++jt){
			cout << jt->first << " ";
		} cout << " )" << endl;
		}*/
		
	//use this to set up a new Network
	/*map< int, vector<link> > lk; 
	for( auto it = mp.begin(); it != mp.end(); ++it){
		auto jt = it->second.begin();
		lk[ it->first ] = vector<link>{ make_pair(jt->first, jt->second) };
		//cout << "induced_graph:: node " << it->first << " nbrs ( " << jt->first << " ";
		++jt;
		for( ;jt != it->second.end(); ++jt){
			lk[ it->first ].push_back( make_pair(jt->first, jt->second) );
			//cout << jt->first << " ";
		} //cout << ")" << endl;
	}
	
	return WeightedNetwork(lk);*/
	
}


WeightedNetwork::WeightedNetwork(const char* filename){
	network_type="weighted_network";	
	read(filename);
	//compute_indices();
	//reindex();
	
}

WeightedNetwork::WeightedNetwork(adjacency_map &mp) : Network(mp){
	network_type="weighted_network";	
}


WeightedNetwork::WeightedNetwork(const Network &G) : Network(G){
	network_type="weighted_network";	
}




#endif
