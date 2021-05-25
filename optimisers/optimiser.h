#ifndef OPTIMISER_H
#define OPTIMISER_H

#include "../qualities/quality.h"
#include "../qualities/modularity.h"
#include "../qualities/barber_modularity.h"
#include "../qualities/projected_modularity.h"
#include <stdlib.h>

using namespace std; 

//interface for community optimiser
template <typename T>
class Optimiser {
	public:
		
		T &Q; 
		string optimiser_type;
		double val;
		int num_nodes;
		int num_communities;
		double eps;
		map<int, int> labels;
		vector<int> new_to_old;
		
		void setup(T& _Q);
		Optimiser(T& _Q);

		virtual double optimise() = 0;
		void shuffle(vector<int> &x);
		void print_basic(bool original_labels=true);
		void print_labels();
		void fixupQ();
		double greedy_update(); 
		void assign_labels(vector< vector<int> > community_maps);
};

///compute final communtiy assignment after louvain steps
template <typename T> void Optimiser<T>::assign_labels(vector< vector<int> > community_maps){
 	labels.clear();
	for(unsigned i=0; i<new_to_old.size(); ++i){ labels[i] = i; }
	
	//community assignment
	for(unsigned i=0; i<community_maps.size(); ++i){
		for(unsigned j=0; j<labels.size(); ++j){
			labels[ j ] = community_maps[i][ labels[ j ] ];
		}
	}

}

///Q fix inconsistencies in Q caused by node swapping
template <typename T> void Optimiser<T>::fixupQ(){
	Q.setup_comms( Q.node_to_comm );
	Q.init( );
	val = Q.eval();
}

///shuffle the elements of x in place
template <typename T> void Optimiser<T>::shuffle(vector<int> &x){
	
	int n = x.size();
	for(int i=0; i<n; ++i){
		int ridx = rand()%n;
		int tmp = x[i];
		x[i] = x[ridx];
		x[ridx] = tmp;
	}
	
}

template <typename T> void Optimiser<T>::print_basic(bool original_labels){

	cout << "this is a " << optimiser_type << " optimiser" << endl;
	if( labels.size() == 0 ){ cout << "No labels found. Run optimise()?" << endl; }
	else{
		for(auto i=labels.begin(); i != labels.end(); ++i){
			int node = (original_labels) ? new_to_old[i->first] : i->first;
			cout << "node(" << node << ") in community " << i->second << endl;
		}	
		cout << "optimal value " << val << endl;
	}
}
///community assignment output
template <typename T> void Optimiser<T>::print_labels(){

	if( labels.size() == 0 ){ cout << "No labels found. Run optimise()?" << endl; }
	else{
		for(auto i=labels.begin(); i != labels.end(); ++i){
			cout << new_to_old[i->first] << " " << i->second << endl;
		}	
	}
}
template <typename T> void Optimiser<T>::setup(T& _Q){

	num_nodes = Q.G.num_nodes;
	num_communities = Q.num_communities;
	
}

template <typename T> Optimiser<T>::Optimiser(T& _Q) : Q(_Q) {
	
	optimiser_type = "optimiser";
	setup(_Q);
	eps = 1e-8;
	new_to_old = Q.G.new_to_old;
	
}

///Equivalent to one_level of a Louvain typr algorithm
template <typename T> double Optimiser<T>::greedy_update(){
	
	double dq;
	int moves;
	int passes = 0;
	
	//D(cout << optimiser_type << "::" << "GREEDY_UPDATE " << endl);
	
	do {
		moves = 0;
		//random node visit order
		vector<int> node_list(num_nodes); for(int i=0; i<num_nodes; ++i){ node_list[i] = i; } shuffle(node_list);

		vector<int> starting = Q.node_to_comm;
		
		for (int i = 0 ; i < num_nodes; ++i) {
			
			int node = node_list[i];
			int comm = Q.node_to_comm[node];
			int d = Q.G.degrees[node];
			//get links from node to other communities
			vector<double> dc = Q.links_to_comms(node);
			

			//random neighbour visit order
			vector<int> nbr_list(Q.G.get_nbrs[node].size()); for(unsigned i=0; i<nbr_list.size(); ++i){ nbr_list[i] = Q.G.get_nbrs[node][i].first; } shuffle(nbr_list);

			//pop node from comm
			Q.remove(node, comm, dc[comm]);
			int best_comm = comm;
			dq = 0;			
			
			for(unsigned j=0; j<nbr_list.size(); ++j){ 
				int nbr = nbr_list[j];
				if( node == nbr ){ continue; } //self or already in same comunity as nbr, dq=0
				
				double g = Q.gain(node, Q.node_to_comm[nbr], dc[Q.node_to_comm[nbr]]); //calculate quality change from adding 'node' to community of 'nbr'

				if( g > dq ){ //best change seen so far
					best_comm = Q.node_to_comm[nbr];
					dq = g;
				}
			
			}
			//push node to comm
 			Q.insert(node, best_comm, dc[best_comm]);
 			//check if something changed
			if (best_comm!=comm){++moves;}
			
		} //finish node loop
		

		//D(cout << optimiser_type << "::" << "Done " << moves << " moves" << endl);
		
		//for(unsigned i=0; i<starting.size(); ++i){
		//	cout << optimiser_type << "::" << i << " comm " << starting[i] << " -> " << Q.node_to_comm[i] << endl;
		//}
		
			
		++passes;
		
  //if we moved a node we must have improved the quality
  } while (moves>0);
  
  //D(cout << optimiser_type << "::" << "Done " << passes << " passes" << endl);
  
  //Q comm_to_nodes is now wrong and the in and tot values are also wrong so fix them
  double old_val = val;
  fixupQ();
  //D(cout << optimiser_type << "::" << "Old Q " << old_val << " new Q " << val << endl);
  
  return val-old_val;
	
}


#endif
