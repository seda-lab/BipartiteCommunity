#ifndef MODULARITY_H
#define MODULARITY_H

#include "quality.h"

using namespace std; 

template<typename T>
class Modularity : public Quality<T> {
	public:
	
		//This is necessary for a templated class
		//Doing this way avoids having to sue this-> in the method definitions
		using Quality<T>::num_communities;
		using Quality<T>::val;
		using Quality<T>::gamma;
		using Quality<T>::node_to_comm;
		using Quality<T>::comm_to_node;
		using Quality<T>::weights;
		using Quality<T>::new_to_old;
		using Quality<T>::in;
		using Quality<T>::tot;
		using Quality<T>::G; 
		using Quality<T>::quality_type;
		
		Modularity(T& _G);
		Modularity(T& _G, vector<int> &labels);
	
		void init();
		double eval();
		void remove(int node, int comm, double degree_in_comm);
		void insert(int node, int comm, double degree_in_comm);
		void join(int comm1, int comm2, double v12);

		double gain(int node, int comm, double degree_in_comm);
		double join_gain(int comm1, int comm2, double v12);
};


template <typename T>
Modularity<T>::Modularity(T& _G) : Quality<T>(_G){
	quality_type = "modularity";	
	Modularity<T>::init();
}

template <typename T>
Modularity<T>::Modularity(T& _G, vector<int> &labels) : Quality<T>(_G, labels){
	quality_type = "modularity";	
	Modularity<T>::init();
}

///set up in and tot vectors
template <typename T>
void Modularity<T>::init(){

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

	//loop over nodes
	for(unsigned i=0; i<G.degrees.size(); ++i){
		tot[ node_to_comm[i] ] += G.degrees[i];
	}
	
}

template <typename T>
double Modularity<T>::eval(){
	
	double q = 0;
	for(int i=0; i<num_communities; ++i){
		if(comm_to_node.size()>0){		//tag for empty community 
			//D( cout << "\n" << quality_type << " eval:: community " << i << " in=" << in[i] << " tot=" << tot[i] << " 2E=" << G.num_links << endl );
			q += in[i] - gamma*tot[i]*tot[i]/G.num_links;
		}
	}
	return q/G.num_links;
}

template <typename T>
void Modularity<T>::remove(int node, int comm, double degree_in_comm){

	//D( cout << quality_type << " remove:: node " << node << " from " << comm << " knc= " << degree_in_comm << " loop = " << G.selfloop(node) << " kn = " << G.degrees[node] <<
	//" in[" << comm << "] " << in[comm] << "->" << in[comm] - (2*degree_in_comm + G.selfloop(node)) << "tot[" << comm << "] " << tot[comm] << "->" << tot[comm] - G.degrees[node] << endl );

	in[comm]  -= 2*degree_in_comm + G.selfloop(node);
	tot[comm] -= G.degrees[node];
	node_to_comm[node] = -1; //temporary tag, comm_to_node is now inconsistent

}

template <typename T>
void Modularity<T>::insert(int node, int comm, double degree_in_comm){

	//D( cout << quality_type << " insert:: node " << node << " into " << comm << " knc= " << degree_in_comm << " loop = " << G.selfloop(node) << " kn = " << G.degrees[node] <<
	//" in[" << comm << "] " << in[comm] << "->" << in[comm] + (2*degree_in_comm + G.selfloop(node)) << "tot[" << comm << "] " << tot[comm] << "->" << tot[comm] + G.degrees[node] << endl );

	in[comm]  += 2*degree_in_comm + G.selfloop(node);
	tot[comm] += G.degrees[node];
	node_to_comm[node] = comm; //temporary tag, comm_to_node is now inconsistent
  
}

template <typename T>
double Modularity<T>::gain(int node, int comm, double degree_in_comm) {

	//D( cout << quality_type << " gain:: node " << node << " into " << comm << " knc= " << degree_in_comm << " kn = " << G.degrees[node] <<
	//" degree[" << node << "] = " << G.degrees[node] << " tot[" << comm << "] = " << tot[comm] << " dQ = " << (degree_in_comm - tot[comm]*G.degrees[node]/G.num_links) << endl );

	return (degree_in_comm - gamma*tot[comm]*G.degrees[node]/G.num_links); //*(2/G.num_links)
}


template <typename T>
void Modularity<T>::join(int comm1, int comm2, double v12){

	//D(cout << quality_type << "join " << comm1 << " into " << comm2 << " v12 = " << v12 << endl);

	in[comm2]  += in[comm1] + 2*v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;
	
	Modularity<T>::qjoin(comm1, comm2);
  
}


template <typename T>
double Modularity<T>::join_gain(int comm1, int comm2, double v12) {

	//D( cout << quality_type << " agg_gain:: comm1 " << comm1 << " into " << comm2 << " v12 = " << v12  
	//<< " tot[" << comm1 << "] = " << tot[comm1] << " tot[" << comm2 << "] = " << tot[comm2] <<
	//" dQ = " << (2*v12 - 2*tot[comm1]*tot[comm2]/G.num_links) );

	return (v12 - gamma*tot[comm1]*tot[comm2]/G.num_links); //*(2/G.num_links)*(2)

}
#endif
