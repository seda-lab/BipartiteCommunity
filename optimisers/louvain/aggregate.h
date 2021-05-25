#ifndef AGGREGATE_H
#define AGGREGATE_H

#include <unordered_map>
#include <map>
#include "../optimiser.h"

using namespace std; 

template <typename T>
class Aggregate : public Optimiser<T> {
	public:

		bool max_agg; //join communities when mod change is zero
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		//functions
		using Optimiser<T>::fixupQ;
	
		unordered_map< int, unordered_map<int, double> > R; //R[c][d] = connection strength from community c to community d
		unordered_map< int, unordered_map<int, double> > Gain; //Gain[c][d] = quality gain joining community c to community d
	
		Aggregate(T& Q, bool _ma);
		
		double setup_gain(int &best_comm1, int &best_comm2);
		inline void update_r(int &best_comm1, int &best_comm2);
		inline void update_gain(int &best_comm2);
		double optimise();
};


template <typename T> 
Aggregate<T>::Aggregate(T& _Q, bool _ma) : Optimiser<T>(_Q) {
	optimiser_type = "aggregate";
	max_agg = _ma;
	
	//this won't deal with weird labelling
	WeightedNetwork tmp = Q.G.induced_graph( Q.node_to_comm ); //compute the community to community coupling
	
	for(unsigned i=0; i< tmp.get_nbrs.size(); ++i){
		for(unsigned j=0; j<tmp.get_nbrs[i].size(); ++j){
			//only want inter-community links
			if(int(i)!=tmp.get_nbrs[i][j].first && R[tmp.get_nbrs[i][j].first].find(i) == R[tmp.get_nbrs[i][j].first].end() && R[i].find(tmp.get_nbrs[i][j].first) == R[i].end() ){		
				R[i][ tmp.get_nbrs[i][j].first ] = tmp.get_nbrs[i][j].second;
				R[tmp.get_nbrs[i][j].first][ i ] = tmp.get_nbrs[i][j].second;
			}
		} 
	}
	
	//D(cout << optimiser_type << ":: " << "R=" << endl);
	/*for(auto it = R.begin(); it != R.end(); ++it){
		for(auto j = it->second.begin(); j != it->second.end(); ++j){
			cout << "(" << it->first << "," << j->first << " : " << j->second << ")" << endl;
		}
	}*/

}


template <typename T> 
double Aggregate<T>::setup_gain(int &best_comm1, int &best_comm2){

  double max_gain = 0;

  for(int comm1=0; comm1<num_communities; ++comm1){
	for(auto j=R[comm1].begin(); j!=R[comm1].end(); ++j){
		
		int comm2 = j->first;
		double jg = Q.join_gain(comm1, comm2, j->second);
		Gain[comm1][comm2] = jg;
		//D(cout << optimiser_type << ":: " << "Gain " << comm1 << "->" << comm2 << " " << Gain[comm1][comm2] << endl);  

		if( jg > max_gain ){
			max_gain = jg; 
			best_comm1 = comm1;
			best_comm2 = comm2;
		}

	}  
  } 
  
  //D(cout << optimiser_type << ":: " << "Max Gain " << best_comm1 << "->" << best_comm2 << " " << max_gain << endl);  
  return max_gain;
	
}

template <typename T> 
inline void Aggregate<T>::update_r(int &comm1, int &comm2){
	
	R[comm2].erase(comm1);
	Gain[comm2].erase(comm1);
	
	R[comm1].erase(comm2);
	Gain[comm1].erase(comm2);
		
	for(auto it=R[comm1].begin(); it != R[comm1].end(); ++it){ //nother eighbours of comm1 are now neighbours of comm2
		int c = it->first;	

		if( R[comm2].find(c) == R[comm2].end() ){ //comm2 was not previously linked to c
			R[comm2][c] = it->second;
			R[c][comm2] = it->second;
		} else {								  //comm2 was previously linked with c	
			R[comm2][c] += it->second;
			R[c][comm2] += it->second;				
		}	

		R[c].erase(comm1);
		Gain[c].erase(comm1); //don't really need to do this

	}

	R.erase(comm1); //comm1 no longer exists
	Gain.erase(comm1); 
	--num_communities;

	/*D(cout << optimiser_type << ":: " << "R=" << endl);
	for(auto it = R.begin(); it != R.end(); ++it){
		for(auto j = it->second.begin(); j != it->second.end(); ++j){
			cout << "(" << it->first << "," << j->first << " : " << j->second << ")" << endl;
		}
	}*/
	  
}

template <typename T> 
inline void Aggregate<T>::update_gain(int &best_comm2){ 

	for(auto it=R[best_comm2].begin(); it != R[best_comm2].end(); ++it){
		int c = it->first;
		if(best_comm2 == c){ continue; }
		double jg = Q.join_gain(best_comm2, c, it->second);
		Gain[best_comm2][c] = jg;
		Gain[c][best_comm2] = jg;
		//D(cout << optimiser_type << ":: " << "Gain " << c << "->" << best_comm2 << "=" << Gain[c][best_comm2] << endl);
	}  

}

template <typename T> 
double Aggregate<T>::optimise() {

	int joins = 0;
	int best_comm1;
	int best_comm2;
	double max_gain = setup_gain(best_comm1, best_comm2);
	
	while( (max_agg) ? (max_gain >= 0) : (max_gain > eps)  ){
		
		//update community labels
		//D(cout << optimiser_type << ":: " << "JOIN " << best_comm1 << "->" << best_comm2 << endl);		
		Q.join(best_comm1, best_comm2, R[best_comm1][best_comm2]);//join comm1 and comm2, we get rid of comm1 and keep comm2!
		++joins;

		//D(cout << optimiser_type << ":: " << "UPDATE R" << endl);			  
		update_r(best_comm1, best_comm2);
		//D(cout << optimiser_type << ":: " << "UPDATE GAIN" << endl);			  		
		update_gain(best_comm2);
		
		//D(cout << optimiser_type << ":: " << "FIND MAX GAIN" << endl);			  				
		max_gain = 0;
		best_comm1 = -1;
		best_comm2 = -1;
		for(auto i=R.begin(); i!=R.end(); ++i){
			int comm1 = i->first;
			for(auto j=i->second.begin(); j != i->second.end(); ++j){
				int comm2 = j->first;
				if( comm1 < comm2 ){ //symmetry
					double jg = Gain[comm1][comm2];
					if( (max_agg) ? (jg >= max_gain) : (jg > max_gain) ){
						max_gain = jg; 
						best_comm1 = comm1;
						best_comm2 = comm2;
					}
				}
			}		  
		}
		if( best_comm1 == best_comm2 ){ break; }
		
		//D(cout << optimiser_type << ":: " << joins << " Max gain = " << max_gain << " for " << best_comm1 << "->" << best_comm2 << endl);	

	}
	//cout << "Aggregation complete after " << joins << " joins" << endl;
	fixupQ();
	//D(cout << optimiser_type << ":: " << "Q = " << endl);	
	//Q.print_basic(false);
	for(unsigned i=0; i<Q.node_to_comm.size(); ++i){labels[i] = Q.node_to_comm[i];}
	
	return val;
}


#endif
