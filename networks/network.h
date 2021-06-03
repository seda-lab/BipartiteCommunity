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
#include <unordered_map>
#include <unordered_set>

using namespace std; 

typedef pair<int, double> link;
typedef unordered_map< int, unordered_map<int, double> > adjacency_map;

///Interface for network objects
class Network {
	public:
		int num_nodes; 		
		double num_links;

		unordered_map< int, double > degrees;			
		//vector< vector<link> > get_nbrs; 	///neighbour list
		adjacency_map get_nbrs;
		
		//depreciated
		//vector< double > cum_degrees;  		///cumulative degree distribution cum_degree[i] = degree[0] + degree[1] + ... + degree[i]; 
		//vector< link > neighbours; 			///edge list
		
		vector<int> new_to_old; 			///mapping between new and old node ids
		//depreciated
		map<int, int> old_to_new;			///mapping between old and new node ids
		
		string network_type;				///Identifier
		
		Network(){
			num_nodes = 0;
			num_links = 0;
			//degrees.resize(0);
			//cum_degrees.resize(0);
			//neighbours.resize(0);
			//get_nbrs.resize(0);
			degrees.clear();
			get_nbrs.clear();
			network_type="network";
		}		
		
		//Network(map< int, vector<link> > &lk);
		Network(adjacency_map &mp);
		void set(adjacency_map &mp);
		void read(const char* filename);

		void compute_degrees();
		void compute_indices();
		void reindex();
		adjacency_map to_map();
		
		void print_basic(bool use_original_ids=true);
		double selfloop(int node);
		Network induced_graph(unordered_map<int, int> &label);
		void join(int i, int j);
		
};


//A_i',a = A_ja + A_ia
//combine node j into i
void Network::join(int i, int j){ 
	if(i != j){
		for(auto &k : get_nbrs[j]){ //k is neighbour of j
			int a = k.first;
			if(a == i){ //connections between i and j become self loops
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){ //does i already have a self loop?
					get_nbrs[i][i] = 2*k.second;
				} else {
					get_nbrs[i][i] += 2*k.second;
				}
				get_nbrs[a].erase(j);  
			} else if(a == j) { //self loops of j remain self loops
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){
					get_nbrs[i][i] = k.second;
				} else {
					get_nbrs[i][i] += k.second;
				}
			} else {
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){
					get_nbrs[i][a] = k.second;
					get_nbrs[a][i] = k.second;
				} else {
					get_nbrs[i][a] += k.second;
					get_nbrs[a][i] += k.second;
				}
				get_nbrs[a].erase(j);  
			}

		}
		degrees[i] += degrees[j];

		get_nbrs.erase(j);  
		degrees.erase(j);
		--num_nodes;
	}
}

Network::Network(adjacency_map &mp){
	network_type="network";		
	set(mp);
}

///Set up a network object from a map
void Network::set(adjacency_map &mp){
	
	get_nbrs = mp;
	num_nodes = mp.size();
	num_links = 0; 
	for(auto it = mp.begin(); it != mp.end(); ++it){ 
		for (auto &j : it->second) {
			num_links += j.second ;
		} 
	}
	//compute_indices();
	//reindex();
	compute_degrees();
}

///return a map representation of the Network
adjacency_map Network::to_map(){
	
	adjacency_map mp;
	
	for(unsigned i=0; i<mp.size(); ++i){ 		
		for( auto &j : get_nbrs[ i ] ){ mp[ new_to_old[i] ][ new_to_old[j.first] ] = j.second; }
	} 	
	
	return mp;
	
}

void Network::compute_degrees(){

	degrees.clear();
	for(auto &i: get_nbrs){
		double d = 0;
		for( auto &j : i.second ){ d += j.second; }
		degrees[i.first] = d;
	}
	
}

///Ensure node labels go from 0 to num_nodes-1. Keep a map to the input labels
void Network::compute_indices(){

	new_to_old.resize(0);
	old_to_new.clear();
	for(auto it = get_nbrs.begin(); it != get_nbrs.end(); ++it){
		old_to_new[it->first] = new_to_old.size();
		new_to_old.push_back(it->first);
	}

}

///compute the node degrees & neighbour list
void Network::reindex(){
			
	//neighbours.resize(0);
	//degrees.resize(0);
	//get_nbrs.resize(num_nodes);
	degrees.clear();
	adjacency_map tmp;
	
	for(unsigned i=0; i<get_nbrs.size(); ++i){ 

		double d = 0;
		for( auto &j : get_nbrs[ new_to_old[i] ] ){ d += j.second; }
		degrees[i] = d;
		//degrees.push_back( d ); 
		//(i==0) ? cum_degrees.push_back( d ) : cum_degrees.push_back( cum_degrees.back() + d ); 
		
		//get_nbrs[i].resize(0);
		for( auto &j : get_nbrs[ new_to_old[i] ] ){
			//neighbours.push_back( make_pair(  old_to_new[ j.first ], j.second )  );
			//get_nbrs[i].push_back( neighbours.back() );
			tmp[i][ old_to_new[ j.first ] ] = j.second;
		} 
		
	}
	
	get_nbrs = tmp; 	
		
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
void Network::read(const char* filename){
	
	ifstream infile(filename);
	if (infile.is_open() != true) {
		cerr << "The file " << filename << " does not exist" << endl;
		exit(1);
	}
  	
	string line;


	num_links = 0;	
	//map< int, vector<link> > lk;
	
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
		get_nbrs[a][b] = w;
		if(a != b){ get_nbrs[b][a] = w; }
		else{ get_nbrs[a][a] += w; }
		/*if( lk.find(a) == lk.end() ){
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
		}*/
		num_links+=2*w;
			
	}
	infile.close();
	//num_nodes = lk.size();
	num_nodes = get_nbrs.size();
	//return lk;	
  	compute_degrees();

  	
}




///return the self edge weight at n. 
double Network::selfloop(int n){
	/*for (auto &j : get_nbrs[n]) {
		if(n == j.first){ return j.second; }
	}*/
	auto it = get_nbrs[n].find(n);
	if(it != get_nbrs[n].end()){
		return it->second;
	}
	return 0;
}

///Use a node labelling to compute an induced network	
Network Network::induced_graph(unordered_map<int, int> &label){
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
	//for(unsigned i=0; i<get_nbrs.size(); ++i){ 
	for(auto &ii : get_nbrs){ 
		int a = ii.first;
		//for( auto &j : get_nbrs[i] ){
		for( auto &j : ii.second ){
			//int a = (use_original_ids) ? new_to_old[i] : i;
			//int b = (use_original_ids) ? new_to_old[ neighbours[cid].first ] : neighbours[cid].first;
			//cout << a << "---[" << neighbours[cid].second << "]---" << b << endl;

			//int b = (use_original_ids) ? new_to_old[ j.first ] : j.first;
			int b = j.first;
			cout << a << "---[" << j.second << "]---" << b << endl;
			
			//++cid;
		}
	}
	//for(unsigned i=0; i<get_nbrs.size(); ++i){ 
	for(auto &ii : get_nbrs){ 
		//int a = (use_original_ids) ? new_to_old[i] : i;
		int a = ii.first;
		cout << "neighbours of " << a << " ( ";
		//for( auto &j : get_nbrs[i] ){		
		for( auto &j : ii.second ){
			int b = j.first; //(use_original_ids) ? new_to_old[ j.first ] : j.first;
			cout << b << " ";
		} cout << ")" << endl;
	}
	double cum_degree = 0;
	//for(unsigned i=0; i<degrees.size(); ++i){ 
	for(auto &i: degrees){ 
		int a = i.first; //(use_original_ids) ? new_to_old[i] : i;
		cum_degree += i.second; //degrees[i];
		cout << "degree " << a << ": " << i.second << " " << cum_degree << endl; 
	}
}


#endif
