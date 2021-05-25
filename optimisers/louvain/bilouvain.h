#ifndef BILOUVAIN_H
#define BILOUVAIN_H


#include "../optimiser.h"
#include "aggregate.h"

using namespace std; 

template<typename T>      
class BiLouvain : public Optimiser<T> {
	public:
		bool aggregate;
		bool max_agg;
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
				

		BiLouvain(T& Q);	
		BiLouvain(T& Q, bool _agg);	
			
		double optimise(); //override
};


template <typename T>
BiLouvain<T>::BiLouvain(T& _Q) : Optimiser<T>(_Q) {
		optimiser_type = "bilouvain";
		aggregate = true;
}
template <typename T>
BiLouvain<T>::BiLouvain(T& _Q, bool _agg) : Optimiser<T>(_Q) {
		optimiser_type = "bilouvain";
		aggregate = _agg;
		max_agg = true;
}

template <typename T>
double BiLouvain<T>::optimise() {
	
	double dq = 0;
	int level = 0;
	val = Q.eval();
	vector< vector<int> > community_maps;
	bool left = (2*Q.G.B > Q.G.num_nodes);
	BipartiteNetwork Original( Q.G );

	//D(cout << optimiser_type << "::" << "Starting at Q = " << val << endl);
	
	do {

		//find new community labels	
		dq = greedy_update(); 
		
		//save community assignment
		community_maps.push_back( Q.node_to_comm );

		//use communities to make induced graph
		BipartiteNetwork BI = Q.G.induced_graph_side( Q.node_to_comm , left ); //condense big size
		//D(cout << optimiser_type << "::" << "Induced Graph = " << endl);
		//BI.print_basic(false);		

		//reset quality function with GI
		Q.setup_network(BI);
		setup(Q);		
		left = !left;
		//D(cout << optimiser_type << "::" << "Modularity = " << endl);		
		//Q.print_basic();
		
		++level;
		//cout << "level " << level << " dq = " << dq << endl;
		
	} while(dq > eps);
	
	assign_labels(community_maps);

	
	if(aggregate){
		//D(cout << optimiser_type << "::" << "Do aggregate = " << endl);		
		//Original.print_basic();
		Q.setup_network(Original);
		vector<int> tlabels( labels.size() );  for(auto &i : labels){ tlabels[ i.first ] = i.second; }
		Q.setup_comms(tlabels);
		Q.init();
		setup(Q);		
		//Q.print_basic(false);
		//D(cout << optimiser_type << "::" << "CALLING AGG................." << endl);		
		Aggregate<T> agg( Q , max_agg );
		agg.optimise();
		//agg.print_basic();
		//agg.print_labels();	
		labels.clear();
		for(auto &i: agg.labels){ labels[ i.first ] = i.second; }
		val = agg.val;
	}
	
	
	return val;
}

#endif
