#ifndef PROJECTION_H
#define PROJECTION_H

#include "network.h"
#include "weighted_network.h"
#include "bipartite_network.h"
#include <algorithm>

using namespace std; 


class Projector {
	public:
		
		bool include_loops;	///should we include self loops (nightmare)
		bool left;			///project on left?
		
		int start;			///Start id of nodes to be kept
		int end;			///end id of nodes to be kept

		int alt_start;		///Start id of nodes to be projected out
		int alt_end;		///end id of nodes to be projected out
				
		int inc;			///0 or 1 if we include or exclude loops
		BipartiteNetwork &G;	///The network we are projecting
		WeightedNetwork GP;		///The projection
		
		
		///Figure out start and end indices.
		void setup_counter(){ 
			if(left){
				start = 0; end = G.B;
				alt_start = G.B; alt_end = G.num_nodes; 
			} else {
				alt_start = 0; alt_end = G.B;
				start = G.B; end = G.num_nodes; 
			}
			if(include_loops){
				inc = 0;
			} else {
				inc = 1;
			}
		}
		Projector(BipartiteNetwork &_G) : G(_G) { left = true; include_loops = false; setup_counter(); }
		Projector(BipartiteNetwork &_G, bool _left, bool _include_loops) : G(_G) { left = _left; include_loops = _include_loops; setup_counter(); }
		
		
		///Compute B.B^T
		map<int, vector<link> > compute_common(vector<int> &weight){
			
			map<int, vector<link> > lk;	
			for(int i=start; i<end; ++i){
				for(int j=i+inc; j<end; ++j){
					double w = 0;
					for(auto &k1 : G.get_nbrs[i]){
					for(auto &k2 : G.get_nbrs[j]){
						if(k1.first == k2.first){ 
							w += k1.second * k2.second / weight[k1.first]; 
						}
					}} 

					if( w > 0 ){
						int a = G.new_to_old[i];
						int b = G.new_to_old[j];
						if( lk.find(a) == lk.end() ){
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
						}
					}
					
				}
			}
			
			if(include_loops){

				for(map<int, vector<link> >::iterator it = lk.begin(); it != lk.end(); ++it){ 
					for(auto j = it->second.begin(); j != it->second.end(); ++j){
						if(j->first == it->first){j->second /= 2;}
					}
				}
				
			}
			
			return lk;
			
		}
		
		virtual WeightedNetwork project() = 0;
		virtual double dsum() = 0;
};


///Binary projection
class BinaryProjector : public Projector {
	public:
	
		BinaryProjector(BipartiteNetwork &_G) : Projector(_G) {}
		BinaryProjector(BipartiteNetwork &_G, bool _left, bool _inc) : Projector(_G, _left, _inc) {}
		
		double dsum(){
			double sum = 0;
			for(int i=alt_start; i<alt_end; ++i){ sum += G.degrees[i]*(G.degrees[i]-inc); }
			return sum;		
		}
		WeightedNetwork project(){
			vector<int> weight(G.num_nodes, 1);	
			map<int, vector<link> > lk = compute_common(weight);
			for(map<int, vector<link> >::iterator it = lk.begin(); it != lk.end(); ++it){ 
				for(auto &j : it->second){ 
					j.second = 1;
				} 
			}
			GP.set(lk);
			return GP;
		}

};

///Weighted projection
class WeightedProjector : public Projector {
	public:

		WeightedProjector(BipartiteNetwork &_G) : Projector(_G) {}
		WeightedProjector(BipartiteNetwork &_G, bool _left, bool _inc) : Projector(_G, _left, _inc) {}
	
		double dsum(){
			double sum = 0;
			for(int i=alt_start; i<alt_end; ++i){  sum += G.degrees[i]*(G.degrees[i]-inc); }
			return sum;		
		}	
		WeightedNetwork project(){
			vector<int> weight(G.num_nodes, 1);
			map<int, vector<link> > lk = compute_common(weight);
			GP.set(lk);
			return GP;
		}
};

///Newman projection
class HyperbolicProjector : public Projector {
	public:

		HyperbolicProjector(BipartiteNetwork &_G) : Projector(_G) {}
		HyperbolicProjector(BipartiteNetwork &_G, bool _left, bool _inc) : Projector(_G, _left, _inc) {}
			
		double dsum(){
			double sum = 0;
			for(int i=alt_start; i<alt_end; ++i){ sum += G.degrees[i];}
			return sum;				
		}		
		WeightedNetwork project(){
			vector<int> weight(G.num_nodes);
			for(unsigned i=0; i<weight.size(); ++i){ weight[i] = G.degrees[i]-inc;   }

			map<int, vector<link> > lk = compute_common(weight);
			GP.set(lk);
			return GP;
		}
};

#endif
