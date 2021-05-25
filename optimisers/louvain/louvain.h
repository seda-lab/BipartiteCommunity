#ifndef LOUVAIN_H
#define LOUVAIN_H

#include "../optimiser.h"

using namespace std; 

template <typename T>
class Louvain : public Optimiser<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		//functions
		using Optimiser<T>::setup;
		using Optimiser<T>::shuffle;
		using Optimiser<T>::print_basic;
		using Optimiser<T>::print_labels;
		using Optimiser<T>::fixupQ;
		using Optimiser<T>::greedy_update; 
		using Optimiser<T>::assign_labels; 
		
		
		Louvain(T& Q);	
		double optimise();
};

template <typename T>
Louvain<T>::Louvain(T& _Q) : Optimiser<T>(_Q) {
		optimiser_type = "louvain";
}

template <typename T>
double Louvain<T>::optimise() {
	
	double dq = 0;
	int level = 0;
	val = Q.eval();
	vector< vector<int> > community_maps;

	//D(cout << optimiser_type << "::" << "Starting at Q = " << val << endl);
	
	do{

		//find new community labels	
		dq = greedy_update(); 
		
		//save community assignment
		community_maps.push_back( Q.node_to_comm );

		//use communities to make induced graph
		WeightedNetwork GI = Q.G.induced_graph( Q.node_to_comm ); 
		//D(cout << optimiser_type << "::" << "Induced Graph = " << endl);
		//GI.print_basic(false);		

		//reset quality function with GI
		Q.setup_network(GI);
		setup(Q);
		//D(cout << optimiser_type << "::" << "Modularity = " << endl);		
		//Q.print_basic();

		++level;
		//D(cout << optimiser_type << "::" << "Level " << level << " dq = " << dq << endl);
		
	} while(dq > eps);
	
	assign_labels(community_maps);
	return val;
}


#endif
