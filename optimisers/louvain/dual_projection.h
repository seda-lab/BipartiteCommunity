#ifndef DUALPROJECTION_H
#define DUALPROJECTION_H


#include "../optimiser.h"
#include "aggregate.h"
#include "projected_louvain.h"

using namespace std; 

template <typename T, typename S>
class DualProjection : public Optimiser<T> {
	public:
		//data
		bool max_agg;
		
		using Optimiser<T>::Q; //the joining Qfn
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		//functions
		using Optimiser<T>::setup;
		using Optimiser<T>::print_basic;
		using Optimiser<T>::print_labels;
		using Optimiser<T>::fixupQ;
		using Optimiser<T>::assign_labels; 
		
		DualProjection(T& _Q);		
		double optimise();
};

template <typename T, typename S>
DualProjection<T, S>::DualProjection(T& _Q) : Optimiser<T>(_Q) {
		optimiser_type = "dual_projection";
		max_agg = true;
}		

template <typename T, typename S>
double DualProjection<T, S>::optimise() {
	
	double dq = 0;
	int level = 0;
	val = Q.eval();
	//BipartiteNetwork Original( Q.G ); //Want an unmodified copy of this

	D(cout << optimiser_type << "::" << "Starting at Q = " << val << endl);
	
	
	//Left communities
	num_communities = 0;
	do{
		BipartiteNetwork B( Q.G );
		S P(B, true, false ); 
		WeightedNetwork W = P.project();
		//W.print_basic();
		ProjectedModularity<S> QP(P, W);
				
		ProjectedLouvain< ProjectedModularity<S> > L(QP);
		L.optimise();
		//L.print_basic();
		for(auto i=L.labels.begin(); i != L.labels.end(); ++i){
			//cout << L.new_to_old[i->first] << " -> " << i->second + num_communities << endl;
			Q.node_to_comm[ L.new_to_old[i->first] ] = i->second + num_communities;
		}	
		num_communities += L.num_communities;
	} while(0);
	
	
	//Right communities
	do{
		BipartiteNetwork B( Q.G );
		S P(B, false, false ); 
		WeightedNetwork W = P.project();
		//W.print_basic();
		ProjectedModularity<S> QP(P, W);
				
		ProjectedLouvain< ProjectedModularity<S> > L(QP);
		L.optimise();
		//L.print_basic();
		for(auto i=L.labels.begin(); i != L.labels.end(); ++i){
			//cout << L.new_to_old[i->first] << " -> " << i->second + num_communities << endl;
			Q.node_to_comm[ L.new_to_old[i->first] ] = i->second + num_communities;
		}		

		num_communities += L.num_communities;
	} while(0);


	//Aggregate
	//D(cout << optimiser_type << "::" << "Do aggregate = " << endl);		
	Q.setup_comms( Q.node_to_comm );
	Q.init();
	setup(Q);		
	//Q.G.print_basic();
	//Q.print_basic();
	//D(cout << optimiser_type << "::" << "CALLING AGG................." << endl);		
	Aggregate<T> agg( Q , max_agg );
	agg.optimise();	
	//agg.print_basic();

	setup(Q);	
	labels.clear();
	for(auto &i: agg.labels){ labels[ i.first ] = i.second; }
	val = agg.val;
	
	return val;
}


#endif
