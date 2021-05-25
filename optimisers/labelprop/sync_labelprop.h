#ifndef SYNCLABELPROP_H
#define SYNCLABELPROP_H

#include <unordered_set>
#include <set>
#include <algorithm>

#include "labelprop.h"

using namespace std; 

template <typename T>
class SyncLabelProp : public LabelProp<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		using LabelProp<T>::high_freq_labels; 
		//functions
		using Optimiser<T>::fixupQ;		
		using LabelProp<T>::calc_high_freq_labels; 
		using LabelProp<T>::use_ew; 

		
		vector< unordered_set<int> > colouring;
		
		SyncLabelProp(T& Q, bool _use_ew=false);
		
		void colour();
		double optimise();
};

template <typename T>
SyncLabelProp<T>::SyncLabelProp(T& _Q, bool _use_ew) : LabelProp<T>(_Q, _use_ew) {
	optimiser_type = "sync_labelprop";
}

///Synchronous label propagation. All labels updated at same time. 
///This is the same algorithm networkx uses
template <typename T>
double SyncLabelProp<T>::optimise() {
	colour();
	/*for(unsigned i=0; i<colouring.size(); ++i){
		cout << "colour " << i << " (";
		for(auto j=colouring[i].begin(); j!=colouring[i].end(); ++j){
			cout << *j << " ";
		} cout << " )" << endl;
	}*/
	
	bool stable = true;
	int passes = 0;
	
	do{
		
		//labels vector contains colours of all the nodes.
		for(unsigned i=0; i<colouring.size(); ++i){
			//cout << "colour " << i << endl;
			for(auto j=colouring[i].begin(); j!=colouring[i].end(); ++j){
				int node = *j;
				
				//cout << "node" << node << endl;
				pair<int, int> high = calc_high_freq_labels(node);
				
				//Use these to give a new label to node
				if( !high.second ){ 
					//cout << "assign " << node << " to " << high.first << endl;
					Q.node_to_comm[node] = high.first; 
				} /*else {
					cout << node << " already has a good label = " << Q.node_to_comm[node] << endl;
				}*/
				
				
			}
		}

		for(int i=0; i<num_nodes; ++i){
			
			pair<int, int> high = calc_high_freq_labels(i);
			
			/*cout << i << " has label " <<  Q.node_to_comm[i] << " high_freq ( ";
			for(auto it=high_freq_labels[i].begin(); it!=high_freq_labels[i].end(); ++it){ cout << *it << " "; } 
			cout << ")" << endl;*/
			
			if( !high.second ){
				//cout << "node " << i << " has more popular neighbours" << endl; 
				stable = false;
				break;
			}
		}
		++passes;
		
	} while (!stable);


	//cout << "label prop finished in " << passes << " passes" << endl;
	
	//renumber everything
	fixupQ();
	//Q.print_basic();
	
	//final labelling
	for(unsigned i=0; i<Q.node_to_comm.size(); ++i){labels[ i ] = Q.node_to_comm[i];}
	
	return 0;
}

template <typename T>
void SyncLabelProp<T>::colour() {
	//colour all the nodes a different colour.
	//Simple but probably not minimal
	
	//define an order to go through the nodes (by degree), init all labels to -1
	labels.clear();
	vector<int> order(num_nodes); 
	for(int i=0;i<num_nodes;++i){ order[i] = i; labels[i] = -1; } 
	sort (order.begin(), order.end(), [&] (auto& a, auto& b)->bool{ return (Q.G.degrees[a] > Q.G.degrees[b]); }  );
	//for(int i=0;i<num_nodes;++i){ cout << i << " " << order[i] << " " << Q.G.degrees[ order[i] ] << endl; }
	
	int max_col = -1;
    for(int i=0;i<num_nodes;++i){
		int node = order[i];
		
		unordered_set<int> neighbour_colours;
		//cout << "colour:: node " << node << " neighbour cols = { ";
		for(unsigned j=0; j<Q.G.get_nbrs[node].size(); ++j){ 
			neighbour_colours.insert( labels[ Q.G.get_nbrs[node][j].first ] );
			//cout << Q.G.get_nbrs[node][j].first << ":" <<  	labels[ Q.G.get_nbrs[node][j].first ] << ", ";
		} //cout << "}" << endl;
		
		//find the smallest 'colour' not yet seen among the neighbours
		int c = 0;
		for(; ; ++c){
			if(neighbour_colours.find(c) == neighbour_colours.end()){break;}
		}
		//cout << "smallest colour = " << node << " = " << c << endl;
		//is this the largest colour so fat?
		if( c>max_col ){ max_col=c; }
		//relabel node
		labels[node] = c;
	}
	//for(int i=0;i<num_nodes;++i){cout << "node" << i << " is " << labels[i] << endl;}
	colouring.resize(max_col+1);
	for(unsigned i=0;i<labels.size();++i){
		colouring[ labels[i] ].insert( i );
		Q.node_to_comm[i] = labels[i];
	}


}
       
        

#endif
