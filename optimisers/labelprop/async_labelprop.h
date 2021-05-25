#ifndef ASYNCLABELPROP_H
#define ASYNCLABELPROP_H

#include <unordered_set>
#include <set>
#include <algorithm>

#include "labelprop.h"


using namespace std; 


template <typename T>
class AsyncLabelProp : public LabelProp<T> {
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
		using LabelProp<T>::calc_high_freq_labels; 
		using Optimiser<T>::fixupQ; 
		using Optimiser<T>::shuffle; 
		using LabelProp<T>::use_ew; 
		
		AsyncLabelProp(T& Q, bool _use_ew=true);
		
		double optimise();
};

template <typename T>
AsyncLabelProp<T>::AsyncLabelProp(T& _Q, bool _use_ew) : LabelProp<T>(_Q, _use_ew) {
	optimiser_type = "async_labelprop";
}

///pretty dumb algorithm for label prop
template <typename T>
double AsyncLabelProp<T>::optimise() {
	
	int passes = 0;
	int moves;
	do{
		
	
		vector<int> node_list(num_nodes); for(int i=0; i<num_nodes; ++i){ node_list[i] = i; } shuffle(node_list); //random order
		moves = 0;
		for(int i=0; i<num_nodes; ++i){
			int node = node_list[i];
			//cout << "operating on node " << node << endl;
			pair<int, int> high = calc_high_freq_labels(node);
			//cout << "high = " << high.first << "," << high.second << endl;
			if( !high.second ){ 
				//cout << "moving " << Q.node_to_comm[node] << "->" << high.first << endl;
				Q.node_to_comm[node] = high.first; //should be a random choice but this ought to do!
				++moves; 
			}
		}
        //cout << "PASS OVER" << endl;
		++passes;
		
	} while (moves > 0);


	//cout << "async label prop finished in " << passes << " passes" << endl;
	
	//renumber everything
	fixupQ();
	//Q.print_basic(false);
	
	//final labelling
	for(unsigned i=0; i<Q.node_to_comm.size(); ++i){labels[ i ] = Q.node_to_comm[i];}

	
	return 0;
}


#endif
