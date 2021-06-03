#ifndef QUALITY_H
#define QUALITY_H

#include "../networks/network.h"
#include "../networks/weighted_network.h"
//#include "../networks/bipartite_network.h"
//#include "../networks/projector.h"
#include <unordered_set>
#include <unordered_map>

using namespace std; 

template <typename T>
class Quality {
	public:
		int num_communities;				//number of communities
		double val;							//value of quality fn
		double gamma;						//resolution parameter
			
		unordered_map< int, int > node_to_comm;			//node to community map
		unordered_map< int, unordered_set<int> > comm_to_node;	//community to node map
		//vector<double> weights;				//depreciated
		//vector<int> new_to_old;				//labels depreciated?
	
		//the partial sums used by Louvain
		unordered_map<int, double> in;		
		unordered_map<int, double> tot;
		unordered_map<int, double> dc; //neighbour communities
				
		T& G; //this must be initialised at construction
		
		string quality_type;

		void setup_comms(unordered_map<int, int> &labels);		
		void setup_network(T& _G);		
		Quality(T& _G);
		Quality(T& _G, unordered_map<int, int> &labels);
		//Quality(T& _G, unordered_map<int, int> &labels, vector<double> &w);

		virtual void init() = 0;				
		virtual double eval() = 0;
		virtual void remove(int node, int comm, double degree_in_comm) = 0;  //remove a node from a community
		virtual void insert(int node, int comm, double degree_in_comm) = 0;	//insert a node into a communtiy
		virtual double gain(int node, int comm, double degree_in_comm) = 0;	//(unnormalised) gain from joining an isolated node n to a communtiy
		virtual double join_gain(int comm1, int comm2, double v12) = 0;		//(unnormalised_ gain from joining comm1 and comm2
		virtual void join(int comm1, int comm2, double v12) = 0;			//join comm1 and comm2
		void qjoin(int comm1, int comm2);									//shared joining function
		
		unordered_map<int, double> links_to_comms(int node);						
		
		void print_basic(bool use_original_ids=true);
};

///Join comm1 to comm2, using the label of comm2
template <typename T>
void Quality<T>::qjoin(int comm1, int comm2){
	cerr << "not mplemented" << endl; exit(1);

	//for(int i=0; i<comm_to_node[comm1].size(); ++i){ 
	/*for(auto &i : comm_to_node[comm1]){ 
		//node_to_comm[ comm_to_node[comm1][i] ] = comm2; 
		node_to_comm[ i ] = comm2; 
		comm_to_node[comm2].push_back( i ); 
	}
	comm_to_node[comm1].resize(0);*/
}

///Calculate how many links	go from node n to the different communities
template <typename T>
unordered_map<int, double> Quality<T>::links_to_comms(int n){

	cerr << "not implemented" << endl; exit(1);
	/*for(int i=0; i<num_communities; ++i){ dc[i] = 0; }
	for(unsigned i=0; i<G.get_nbrs[n].size(); ++i){
		if(n != G.get_nbrs[n][i].first){
			dc[  node_to_comm[G.get_nbrs[n][i].first]  ]  +=  G.get_nbrs[n][i].second; 
		}
	}*/
	return dc;
	
}

///output basic info
template <typename T>
void Quality<T>::print_basic(bool use_original_ids){

	cout << quality_type << " is the quality function" << endl;
	cout << num_communities << " communities" << endl;
	/*for(unsigned i=0; i<node_to_comm.size(); ++i){ 
		int c = (use_original_ids) ? new_to_old[ node_to_comm[i] ] : node_to_comm[i]; 
		cout << i << " in " << c << endl; 
	}*/
	for(auto &i : node_to_comm){ cout << i.first << " in " << i.second << endl; }
		
	/*for(unsigned i=0; i<new_to_old.size(); ++i){ 
		int c = (use_original_ids) ? new_to_old[i] : i; 
		cout << "community " << c << " = (";
		for(auto &j : comm_to_node[i]){ cout << " " << j; } cout << " )" << endl;
	}*/
	for(auto &i : comm_to_node){
		cout << "community " << i.first << " = (";
		for(auto &j : i.second){
			cout << j << " ";
		} cout << ")" << endl;
	}
		
	
	/*for(unsigned i=0; i<in.size(); ++i){
		cout << "in[" << i << "] = " << in[i]  << endl;
	}
	for(unsigned i=0; i<tot.size(); ++i){
		cout << "tot[" << i << "] = " << tot[i] << endl;
	}*/
	cout << "Current Quality = " << val << endl;
	
}

/*
template <typename T>
Quality<T>::Quality(T& _G, vector<int> &labels, vector<double> &w) : G(_G) { 	
	quality_type = "quality";
	gamma = 1;
	weights = w;
	setup_comms(labels); 
}*/

template <typename T>
Quality<T>::Quality(T& _G, unordered_map<int, int> &labels) : G(_G) { 
	quality_type = "quality";
	gamma = 1; 
	//weights.resize(0); 
	setup_comms(labels); 
}

///Use a network to initialise Q
template <typename T>
void Quality<T>::setup_network(T& _G){ 

	G = _G;

	/*node_to_comm.resize( G.num_nodes );
	comm_to_node.clear(); //resize( G.num_nodes );
	new_to_old.resize( G.num_nodes );
	for(int i=0; i<G.num_nodes; ++i){
		node_to_comm[i] = i;
		comm_to_node[i] = unordered_set<int>{i};
		new_to_old[i] = i;
	}
	num_communities = comm_to_node.size();
	dc.resize(num_communities);*/
	init();
}

template <typename T>
Quality<T>::Quality(T& _G) : G(_G){ 

	quality_type = "quality";
	//weights.resize(0); 
	gamma = 1;
	//vector<int> labels(G.num_nodes); for(unsigned i=0; i<labels.size(); ++i){ labels[i] = i; }
	unordered_map<int, int> labels; for(auto &i : G.get_nbrs){ labels[ i.first ] = i.first; }
	setup_comms( labels );
	
}

///Set up community maps from a labelling
template <typename T>
void Quality<T>::setup_comms(unordered_map<int, int> &labels){ 
	
	if(labels.size() != unsigned(G.num_nodes) ){
		cerr << "Trying to use " << labels.size() << " labels with " << G.num_nodes << " node network" << endl;
		exit(1);
	}
	
	node_to_comm = labels;
	for(auto &i : node_to_comm){
		if(comm_to_node.find(i.second) == comm_to_node.end()){
			comm_to_node[i.second] = unordered_set<int>{ i.first };
		} else {
			comm_to_node[i.second].insert( i.first );
		}
	}
	
	num_communities = comm_to_node.size();
	//dc.resize(num_communities);
	
	/*map<int, vector<int> > mp;
	for(unsigned i=0;i<labels.size(); ++i){
		if( mp.find( labels[i] ) == mp.end() ){ mp[ labels[i] ] = vector<int>{int(i)};  }
		else{ mp[ labels[i] ].push_back(i);  }
	}
	
	comm_to_node.clear(); //resize(0);
	new_to_old.resize(0);
	for(auto it=mp.begin(); it != mp.end(); ++it){
		comm_to_node[ new_to_old.size() ] = unordered_set<int>{};
		for(auto &i: it->second){ comm_to_node[ new_to_old.size() ].insert( i ); }

		new_to_old.push_back( it->first );
	}
	
	node_to_comm.resize(labels.size());*/
	/*for(unsigned i=0; i<comm_to_node.size(); ++i){
		for(unsigned j=0; j<comm_to_node[i].size(); ++j){
			node_to_comm[ comm_to_node[i][j] ] = i;
		}
	}*/
	/*for(auto it = comm_to_node.begin(); it!=comm_to_node.end(); ++it){
		for(auto &j : it->second){
			node_to_comm[ j ] = it->first;
		}
	}*/
	


}



#endif
