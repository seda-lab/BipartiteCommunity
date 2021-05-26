#ifndef PROJECTEDMODULARITY_H
#define PROJECTEDMODULARITY_H

#include "quality.h"

using namespace std; 

template <typename T>
class ProjectedModularity : public Quality< WeightedNetwork > {
	public:
	
		T &P;
		double dsum;
		
		ProjectedModularity(T &_p, WeightedNetwork& _G);
		ProjectedModularity(T &_p, WeightedNetwork& _G, vector<int> &labels);
		
			
		void init();
		double eval();
		void remove(int node, int comm, double degree_in_comm);
		void insert(int node, int comm, double degree_in_comm);
		void join(int comm1, int comm2, double v12);

		double gain(int node, int comm, double degree_in_comm);
		double join_gain(int comm1, int comm2, double v12);
};

template <typename T>
ProjectedModularity<T>::ProjectedModularity(T &_p, WeightedNetwork& _G) :  Quality< WeightedNetwork >(_G), P(_p)  {
	quality_type = "projected_modularity";	
	ProjectedModularity<T>::init();
}

template <typename T>
ProjectedModularity<T>::ProjectedModularity(T &_p, WeightedNetwork& _G, vector<int> &labels) : Quality< WeightedNetwork >(_G, labels), P(_p)  {
	quality_type = "projected_modularity";	
	if(labels.size() != (unsigned)G.num_nodes){ cerr << "need one label per node in projected network for " << quality_type << endl; exit(1); }
	ProjectedModularity<T>::init();
}

template <typename T>
void ProjectedModularity<T>::init(){

	in.resize( num_communities ); 
	tot.resize( num_communities ); 
	for(int i=0; i<num_communities; ++i){ in[i] = 0; tot[i] = 0;}

	
	//loop over edges
	for(unsigned i=0;i<G.get_nbrs.size();++i){
		for(unsigned j=0; j<G.get_nbrs[i].size(); ++j){
			if( node_to_comm[i] == node_to_comm[ G.get_nbrs[i][j].first ] ){
				in[ node_to_comm[i] ] += G.get_nbrs[i][j].second;
			} 
		}
	}

	
	//loop over nodes in the correct bipartite set
	for(unsigned i=0; i<node_to_comm.size(); ++i){
		tot[ node_to_comm[i] ] += P.G.degrees[P.start + i]; 
	}	
	//compute the normalisation factor
	dsum = P.dsum()/( (P.G.num_links/2)*( (P.G.num_links/2)-P.inc) );

}


template <typename T>
double ProjectedModularity<T>::eval(){
	
	double q = 0;
	for(int i=0; i<num_communities; ++i){
		if(comm_to_node.size()>0){		//tag for empty community 
			//D( cout << "\n" << quality_type << " eval:: community " << i << " in=" << in[i] << " dsum=" << dsum << " 2E=" << G.num_links << " += " << in[i] - dsum*tot[i]*tot[i] << endl );	
			q += in[i] - gamma*dsum*tot[i]*tot[i];
		}
	}
	return q/G.num_links;

	return 1;
}

template <typename T>
void ProjectedModularity<T>::remove(int node, int comm, double degree_in_comm){

	//D( cout << quality_type << " remove:: node " << node << " from " << comm << " knc= " << degree_in_comm << " loop= " <<  G.selfloop(node) << " qn = " << P.G.degrees[P.start+node] <<
	//" in[" << comm << "] " << in[comm] << " -> " << in[comm] - (2*degree_in_comm + G.selfloop(node)) << " tot[" << comm << "] " 
	//<< tot[comm] << " -> " << tot[comm] - P.G.degrees[P.start+node] << endl );

	in[comm]  -= 2*degree_in_comm + G.selfloop(node);
	tot[comm] -= P.G.degrees[P.start+node];
	node_to_comm[node] = -1; //temporary tag, comm_to_node is now inconsistent

}

template <typename T>
void ProjectedModularity<T>::insert(int node, int comm, double degree_in_comm){

	//D( cout << quality_type << " insert:: node " << node << " into " << comm << " knc= " << degree_in_comm << " qn = " << P.G.degrees[P.start+node] <<
	//" in[" << comm << "] " << in[comm] << " -> " << in[comm] + (2*degree_in_comm + G.selfloop(node)) << " tot[" << comm << "] " 
	//<< tot[comm] << " -> " << tot[comm] + P.G.degrees[P.start+node] << endl );
	
	in[comm]  += 2*degree_in_comm + G.selfloop(node);
	tot[comm] += P.G.degrees[P.start+node];
	node_to_comm[node] = comm; //temporary tag, comm_to_node is now inconsistent

}

template <typename T>
double ProjectedModularity<T>::gain(int node, int comm, double degree_in_comm) {

	//D( cout << quality_type << " gain:: node " << node << " into " << comm << " knc= " << degree_in_comm << " qn = " << P.G.degrees[P.start+node] <<
	//" dsum= " << dsum << " tot[" << comm << "] = " << tot[comm] << " dQ = " << (degree_in_comm - dsum*tot[comm]*P.G.degrees[P.start+node]) << endl );
	
	return (degree_in_comm - gamma*dsum*tot[comm]*P.G.degrees[P.start+node]); //*(2/G.num_links)
	
}

template <typename T>
void ProjectedModularity<T>::join(int comm1, int comm2, double v12){

	//D(cout << quality_type << "join " << comm1 << " into " << comm2 << " v12 = " << v12 << endl);

	in[comm2]  += in[comm1] + 2*v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;
	
	ProjectedModularity<T>::qjoin(comm1, comm2);
  
}



template <typename T>
double ProjectedModularity<T>::join_gain(int comm1, int comm2, double v12) {

	//D( cout << quality_type << " agg_gain:: comm1 " << comm1 << " into " << comm2 << " v12 = " << v12  
	//<< " tot[" << comm1 << "] = " << tot[comm1] << " tot[" << comm2 << "] = " << tot[comm2] <<
	//" dQ = " << (2*v12 - 2*tot[comm1]*tot[comm2]/G.num_links) );
	
	return (v12 - gamma*dsum*tot[comm1]*tot[comm2]); //*(2/G.num_links)

}

#endif
