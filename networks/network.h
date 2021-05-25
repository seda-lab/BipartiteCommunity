#ifndef NETWORK_H
#define NETWORK_H

#ifdef DEBUG 
	#define D(x) (x)
#else 
	#define D(x) do{}while(0)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std; 

typedef pair<int, double> link;

///Interface for network objects
class Network {
	public:
		int num_nodes; 		
		double num_links;

		vector< double > degrees;			
		vector< vector<link> > get_nbrs; 	///neighbour list

		//depreciated
		vector< double > cum_degrees;  		///cumulative degree distribution cum_degree[i] = degree[0] + degree[1] + ... + degree[i]; 
		vector< link > neighbours; 			///edge list
		
		vector<int> new_to_old; 			///mapping between new and old node ids
		//depreciated
		map<int, int> old_to_new;			///mapping between old and new node ids
		
		string network_type;				///Identifier
		
		Network(){
			num_nodes = 0;
			num_links = 0;
			degrees.resize(0);
			cum_degrees.resize(0);
			neighbours.resize(0);
			get_nbrs.resize(0);
			network_type="network";
		}		
		
		Network(map< int, vector<link> > &lk);
		
		void set(map< int, vector<link> > &lk);
		void compute_indices(map< int, vector<link> > &lk);
		void reindex(map< int, vector<link> > &lk);
		map< int, vector<link> > read(const char* filename);
		map< int, vector<link> > to_map();
		
		void print_basic(bool use_original_ids=true);
		double selfloop(int node);
		Network induced_graph(vector<int> &label);
		
};

Network::Network(map< int, vector<link> > &lk){
	network_type="network";		
	set(lk);
}

///Set up a network object from a map
void Network::set(map< int, vector<link> > &lk){
	num_nodes = lk.size();
	num_links = 0; 
	for(map<int, vector<link> >::iterator it = lk.begin(); it != lk.end(); ++it){ 
		for (auto &j : it->second) {
			num_links += j.second ;
		} 
	}
	compute_indices(lk);
	reindex(lk);
}

///return a map representation of the Network
map< int, vector<link> > Network::to_map(){
	
	map< int, vector<link> > lk;
	
	for(int i=0; i<num_nodes; ++i){ 		
		for( auto &j : get_nbrs[ i ] ){ lk[ new_to_old[i] ].push_back(j); }
	} 	
	
	return lk;
	
}

///Ensure node labels go from 0 to num_nodes-1. Keep a map to the input labels
void Network::compute_indices(map< int, vector<link> > &lk){

	new_to_old.resize(0);
	old_to_new.clear();
	for(map<int, vector<link> >::iterator it = lk.begin(); it != lk.end(); ++it){
		old_to_new[it->first] = new_to_old.size();
		new_to_old.push_back(it->first);
	}

}

///compute the node degrees & neighbour list
void Network::reindex(map< int, vector<link> > &lk){
			
	neighbours.resize(0);
	degrees.resize(0);
	get_nbrs.resize(num_nodes);
	
	for(int i=0; i<num_nodes; ++i){ 


		double d = 0;
		for( auto &j : lk[ new_to_old[i] ] ){ d += j.second; }
		degrees.push_back( d ); 
		(i==0) ? cum_degrees.push_back( d ) : cum_degrees.push_back( cum_degrees.back() + d ); 
		
		get_nbrs[i].resize(0);
		for( auto &j : lk[ new_to_old[i] ] ){
			neighbours.push_back( make_pair(  old_to_new[ j.first ], j.second )  );
			get_nbrs[i].push_back( neighbours.back() );
		} 
		
	} 	
		
}
	

/**
 * Read a network in from a file
 * Expected format
 * weighted = false
 * 0 1
 * 0 2
 * 1 2 ...
 * weighted = true
 * 0 1 1
 * 0 2 1
 * 1 2 3 ...
 * Edges are assumed to be bidirectional!
 * **/
map< int, vector<link> > Network::read(const char* filename){
	
	ifstream infile(filename);
	if (infile.is_open() != true) {
		cerr << "The file " << filename << " does not exist" << endl;
		exit(1);
	}
  	
	string line;


	num_links = 0;	
	map< int, vector<link> > lk;
	
	while(getline(infile, line)) {

		//skip comments and blank lines: %, #, // allowed
		if(line.empty() || (line.find("%") == 0) || (line.find("#") == 0) || (line.find("//") == 0) ) {continue;} 
		
		//read lines
		istringstream iss(line);
		string token;
		int num_tokens = 0;

		int a, b;
		vector<int> n(2);
		double w=1;
	
		while(getline(iss, token, ' ')){
			if(num_tokens < 2){
				n[num_tokens] = atoi(token.c_str());
			} else {
				w = atof(token.c_str());
			}
			++num_tokens;
		}
		
		a = n[0]; b = n[1];
		if( lk.find(a) == lk.end() ){
			lk[a] = vector<link>{ make_pair(b, w) };
		} else {
			lk[a].push_back( make_pair(b, w) ); 						
		}
			
		if( lk.find(b) == lk.end() ){
			lk[b] = vector<link>{ make_pair(a, w) };
		} else {
			if(a != b){
				lk[b].push_back( make_pair(a, w)  ); 		
			} else { 	//only want self loops on neighbour list once, but should contribute 2 to degree!
				auto it = find (lk[b].begin(), lk[b].end(), make_pair(a, w));
				it->second += w;
			}
		}
		num_links+=2*w;
			
	}
	infile.close();
	num_nodes = lk.size();
	return lk;	
  	
  	
}




///return the self edge weight at n. 
double Network::selfloop(int n){
	for (auto &j : get_nbrs[n]) {
		if(n == j.first){ return j.second; }
	}
	return 0;
}

///Use a node labelling to compute an induced network	
Network Network::induced_graph(vector<int> &label){
	cerr << "induced_graph should not be called from Base class Network" << endl;
	exit(1);
	return Network();
}

///print basic infor about the network			
void Network::print_basic(bool use_original_ids){
	cout << "This is a " << network_type << endl;
	cout << num_nodes << " nodes" << endl;
	cout << num_links << " total edge weight" << endl;
	//int cid = 0;
	for(int i=0; i<num_nodes; ++i){ 
		for( auto &j : get_nbrs[i] ){
			int a = (use_original_ids) ? new_to_old[i] : i;
			//int b = (use_original_ids) ? new_to_old[ neighbours[cid].first ] : neighbours[cid].first;
			//cout << a << "---[" << neighbours[cid].second << "]---" << b << endl;

			int b = (use_original_ids) ? new_to_old[ j.first ] : j.first;
			cout << a << "---[" << j.second << "]---" << b << endl;
			
			//++cid;
		}
	}
	for(int i=0; i<num_nodes; ++i){ 
		int a = (use_original_ids) ? new_to_old[i] : i;
		
		cout << "neighbours of " << a << " ( ";
		for( auto &j : get_nbrs[i] ){		
			int b = (use_original_ids) ? new_to_old[ j.first ] : j.first;
			cout << b << " ";
		} cout << ")" << endl;
	}
	for(int i=0; i<num_nodes; ++i){ 
		int a = (use_original_ids) ? new_to_old[i] : i;
		cout << "degree " << a << ": " << degrees[i] << " " << cum_degrees[i] << endl; 
	}
}


#endif
