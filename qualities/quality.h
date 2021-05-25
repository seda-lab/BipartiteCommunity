#ifndef QUALITY_H
#define QUALITY_H

#include "../networks/network.h"
#include "../networks/weighted_network.h"
#include "../networks/bipartite_network.h"
#include "../networks/projector.h"

using namespace std; 

template <typename T>
class Quality {
	public:
		int num_communities;				//number of communities
		double val;							//value of quality fn
			
		vector<int> node_to_comm;			//node to community map
		vector< vector<int> > comm_to_node;	//community to node map
		vector<double> weights;				//depreciated
		vector<int> new_to_old;				//labels depreciated?
	
		//the partial sums used by Louvain
		vector<double> in;		
		vector<double> tot;
		
		T& G; //this must be initialised at construction
		
		string quality_type;

		void setup_comms(vector<int> &labels);		
		void setup_network(T& _G);		
		Quality(T& _G);
		Quality(T& _G, vector<int> &labels);
		Quality(T& _G, vector<int> &labels, vector<double> &w);

		virtual void init() = 0;				
		virtual double eval() = 0;
		virtual void remove(int node, int comm, double degree_in_comm) = 0;  //remove a node from a community
		virtual void insert(int node, int comm, double degree_in_comm) = 0;	//insert a node into a communtiy
		virtual double gain(int node, int comm, double degree_in_comm) = 0;	//(unnormalised) gain from joining an isolated node n to a communtiy
		virtual double join_gain(int comm1, int comm2, double v12) = 0;		//(unnormalised_ gain from joining comm1 and comm2
		virtual void join(int comm1, int comm2, double v12) = 0;			//join comm1 and comm2
		void qjoin(int comm1, int comm2);									//shared joining function
		
		vector<double> links_to_comms(int node);						
		
		void print_basic(bool use_original_ids=true);
};

///Join comm1 to comm2, using the label of comm2
template <typename T>
void Quality<T>::qjoin(int comm1, int comm2){
	

	//for(int i=0; i<comm_to_node[comm1].size(); ++i){ 
	for(auto &i : comm_to_node[comm1]){ 
		//node_to_comm[ comm_to_node[comm1][i] ] = comm2; 
		node_to_comm[ i ] = comm2; 
		comm_to_node[comm2].push_back( i ); 
	}
	comm_to_node[comm1].resize(0);
}

///Calculate how many links	go from node n to the different communities
template <typename T>
vector<double> Quality<T>::links_to_comms(int n){

	vector<double> d(num_communities, 0);
	for(unsigned i=0; i<G.get_nbrs[n].size(); ++i){
		if(n != G.get_nbrs[n][i].first){
			d[  node_to_comm[G.get_nbrs[n][i].first]  ]  +=  G.get_nbrs[n][i].second; 
		}
	}
	return d;
	
}

///output basic info
template <typename T>
void Quality<T>::print_basic(bool use_original_ids){

	cout << quality_type << " is the quality function" << endl;
	cout << num_communities << " communities" << endl;
	for(unsigned i=0; i<node_to_comm.size(); ++i){ 
		int c = (use_original_ids) ? new_to_old[ node_to_comm[i] ] : node_to_comm[i]; 
		cout << i << " in " << c << endl; 
	}
		
	for(unsigned i=0; i<new_to_old.size(); ++i){ 
		int c = (use_original_ids) ? new_to_old[i] : i; 
		cout << "community " << c << " = (";
		for(auto &j : comm_to_node[i]){ cout << " " << j; } cout << " )" << endl;
	}
	
	for(unsigned i=0; i<in.size(); ++i){
		cout << "in[" << i << "] = " << in[i]  << endl;
	}
	for(unsigned i=0; i<tot.size(); ++i){
		cout << "tot[" << i << "] = " << tot[i] << endl;
	}
	cout << "Current Quality = " << eval() << endl;
	
}

template <typename T>
Quality<T>::Quality(T& _G, vector<int> &labels, vector<double> &w) : G(_G) { 	
	quality_type = "quality";
	weights = w;
	setup_comms(labels); 

}

template <typename T>
Quality<T>::Quality(T& _G, vector<int> &labels) : G(_G) { 
	quality_type = "quality";
	weights.resize(0); 
	setup_comms(labels); 

}

///Use a network to initialise Q
template <typename T>
void Quality<T>::setup_network(T& _G){ 

	G = _G;

	node_to_comm.resize( G.num_nodes );
	comm_to_node.resize( G.num_nodes );
	new_to_old.resize( G.num_nodes );
	for(int i=0; i<G.num_nodes; ++i){
		node_to_comm[i] = i;
		comm_to_node[i] = vector<int>{i};
		new_to_old[i] = i;
	}
	num_communities = comm_to_node.size();
	init();
}

template <typename T>
Quality<T>::Quality(T& _G) : G(_G){ 

	quality_type = "quality";
	weights.resize(0); 
	vector<int> labels(G.num_nodes); for(unsigned i=0; i<labels.size(); ++i){ labels[i] = i; }
	setup_comms( labels );
	
}

///Set up community maps from a labelling
template <typename T>
void Quality<T>::setup_comms(vector<int> &labels){ 
	
	if(labels.size() != unsigned(G.num_nodes) ){
		cerr << "Trying to use " << labels.size() << " labels with " << G.num_nodes << " node network" << endl;
		exit(1);
	}
			
	map<int, vector<int> > mp;
	for(unsigned i=0;i<labels.size(); ++i){
		if( mp.find( labels[i] ) == mp.end() ){ mp[ labels[i] ] = vector<int>{int(i)};  }
		else{ mp[ labels[i] ].push_back(i); }
	}
	
	comm_to_node.resize(0);
	new_to_old.resize(0);
	for(auto it=mp.begin(); it != mp.end(); ++it){
		new_to_old.push_back( it->first );
		comm_to_node.push_back( it->second );
	}

	node_to_comm.resize(labels.size());
	for(unsigned i=0; i<comm_to_node.size(); ++i){
		for(unsigned j=0; j<comm_to_node[i].size(); ++j){
			node_to_comm[ comm_to_node[i][j] ] = i;
		}
	}
	num_communities = comm_to_node.size();

}



#endif
