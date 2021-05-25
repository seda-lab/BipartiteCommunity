#ifndef LABELPROP_H
#define LABELPROP_H

#include "../optimiser.h"
#include <unordered_set>
#include <set>
#include <algorithm>

using namespace std; 

template<typename T>      
class LabelProp : public Optimiser<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		
		bool use_ew;
		vector< unordered_set<int> > high_freq_labels; ///for each node what label(s) appear most often among its neighbours?
		
		LabelProp(T& Q, bool _use_ew=false);
		pair<int, int> calc_high_freq_labels(int node);

};

template <typename T>
LabelProp<T>::LabelProp(T& _Q, bool _use_ew) : Optimiser<T>(_Q) {
	optimiser_type = "labelprop";
	use_ew = _use_ew;
	high_freq_labels.resize( num_nodes );	
}      

///get the most frequent labels of node's neighbours
template <typename T>      
pair<int, int> LabelProp<T>::calc_high_freq_labels(int node) {
	
	
	map< int, double > freqs;
	int best_val = -1;
	for(unsigned j=0; j<Q.G.get_nbrs[node].size(); ++j){ 
		int nbr = Q.G.get_nbrs[node][j].first;
		int nbr_label = Q.node_to_comm[nbr];
		double w = (use_ew) ? Q.G.get_nbrs[node][j].second : 1;
		if( freqs.find(nbr_label) == freqs.end()){ freqs[ nbr_label ] = w; }
		else { freqs[ nbr_label ] += w; }
		
		if(freqs[ nbr_label ] > best_val){ best_val = freqs[ nbr_label ]; }
	}

	/*cout << "freqs of " << node << " (";
	for(auto it=freqs.begin(); it!=freqs.end(); ++it){
		cout << it->first << ":" << it->second << " ";
	} cout << ")" << endl;*/
	
	//get all the labels that occur as often as best_val
	unordered_set<int> best;
	int high_label = -1;		//an element of best (node label)
	int current_comm = false;	//is current label one of the ones that occurs most often?

	for(auto it=freqs.begin(); it!=freqs.end(); ++it){
		if( it->second == best_val ){ 
			best.insert(it->first); 
			if( it->first > high_label ){ high_label = it->first; }
			if( it->first == Q.node_to_comm[node] ){ current_comm = true; }
		}
	}
	
	//cout << node << " common neighbour labels ( "; for(auto it=best.begin(); it!=best.end(); ++it){ cout << *it << " ";} cout << ")" << endl;
	
	high_freq_labels[node] = best;
	return make_pair(high_label, current_comm);
}       
        


#endif
