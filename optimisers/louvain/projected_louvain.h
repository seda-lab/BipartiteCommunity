#ifndef PROJECTEDLOUVAIN_H
#define PROJECTEDLOUVAIN_H


#include "../optimiser.h"

using namespace std; 

template <typename T>
class ProjectedLouvain : public Optimiser<T> {
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
		
		ProjectedLouvain(T& _Q);		
		double optimise();
};

template <typename T>
ProjectedLouvain<T>::ProjectedLouvain(T& _Q) : Optimiser<T>(_Q) {
	
	
		optimiser_type = "projected_louvain";
}		

template <typename T>
double ProjectedLouvain<T>::optimise() {
	
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
		//With include loops=0 this is unnecessary otherwise we lose edges weight as we aggregate

		//Original graph
		WeightedNetwork GI = Q.P.GP.induced_graph( Q.node_to_comm ); 
		//D(cout << optimiser_type << "::" << "Induced graph = " << endl);
		//GI.print_basic(false);	

		
		BipartiteNetwork BI = Q.P.G.induced_graph_side(Q.node_to_comm, Q.P.left);
		//D(cout << optimiser_type << "::" << "Bipartite induced graph = " << val << endl);
		//BI.print_basic(false);
		
		//Update Quality	
		Q.P.GP = GI; //GI is NOT the projection of BI in general
		Q.P.G = BI; 
		Q.P.setup_counter();
		Q.setup_network(GI);
		setup(Q);

		//D(cout << optimiser_type << "::" << "Projected Modularity = " << endl);		
		//Q.print_basic();
		
		
		++level;
		//D(cout << optimiser_type << "::" << "Level " << level << " dq = " << dq << endl);
		
	} while(dq > eps);
	
	assign_labels(community_maps);
	return val;
}


#endif
