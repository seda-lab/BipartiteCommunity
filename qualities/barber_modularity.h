#ifndef BARBERMODULARITY_H
#define BARBERMODULARITY_H

#include "quality.h"

using namespace std; 

class BarberModularity : public Quality< BipartiteNetwork > {
	public:
	
		
		
		BarberModularity(BipartiteNetwork& _G);
		BarberModularity(BipartiteNetwork& _G, vector<int> &labels);
	
		void init();
		double eval();
		void remove(int node, int comm, double degree_in_comm);
		void insert(int node, int comm, double degree_in_comm);
		void join(int comm1, int comm2, double v12);

		double gain(int node, int comm, double degree_in_comm);
		double join_gain(int comm1, int comm2, double v12);
};


BarberModularity::BarberModularity(BipartiteNetwork& _G) : Quality<BipartiteNetwork>(_G){
	quality_type = "barber_modularity";	
	BarberModularity::init();
}

BarberModularity::BarberModularity(BipartiteNetwork& _G, vector<int> &labels) : Quality<BipartiteNetwork>(_G, labels){
	quality_type = "barber_modularity";	
	BarberModularity::init();
}

void BarberModularity::init(){

	in.resize( num_communities ); for(int i=0; i<num_communities; ++i){ in[i] = 0; }
	tot.resize( 2*num_communities ); for(int i=0; i<2*num_communities; ++i){ tot[i] = 0; }


	//loop over edges
	for(unsigned i=0;i<G.get_nbrs.size();++i){
		for(unsigned j=0; j<G.get_nbrs[i].size(); ++j){
			if( node_to_comm[i] == node_to_comm[ G.get_nbrs[i][j].first ] ){
				in[ node_to_comm[i] ] += G.get_nbrs[i][j].second;
			} 
		}
	}

	//left degree sums
	for(int i=0; i<G.B; ++i){
		tot[ node_to_comm[i] ] += G.degrees[i];
	}
	//right degree sums	
	for(int i=G.B; i<G.num_nodes; ++i){
		tot[ num_communities + node_to_comm[i] ] += G.degrees[i];
	}

	for(int i=0; i<num_communities; ++i){ in[i] /= 2; } //don't double count in bipartite
	
}


double BarberModularity::eval(){
	
	double q = 0;
	for(int i=0; i<num_communities; ++i){
		if(comm_to_node.size()>0){		//tag for empty community 
			//D( cout << "\n" << quality_type << " eval:: community " << i << " in=" << in[i] << " q=" << tot[i] << " d=" << tot[num_communities+i] << " F=" << G.num_links/2 << 
			//"+=" << in[i] - 2*tot[i]*tot[num_communities+i]/G.num_links <<endl );			
			q += in[i] - 2*tot[i]*tot[num_communities+i]/G.num_links; //factor of 2 is for double counting edges
		}
	}
	return 2*q/G.num_links;
}


void BarberModularity::remove(int node, int comm, double degree_in_comm){


	//D( cout << quality_type << " remove:: node " << node << " from " << comm << " knc= " << degree_in_comm << " kn = " << G.degrees[node] <<
	//" in[" << comm << "] " << in[comm] << " -> " << in[comm] - degree_in_comm << " tot[" << comm + ( (node<G.B)?0:num_communities ) << "] " 
	//<< tot[comm+( (node<G.B)?0:num_communities )] << " -> " << tot[comm+( (node<G.B)?0:num_communities )] - G.degrees[node] << endl );


	in[comm]  -= degree_in_comm;
	if( node < G.B ){
		tot[comm] -= G.degrees[node];
	} else {
		tot[num_communities + comm] -= G.degrees[node];
	}
	
	node_to_comm[node] = -1; //temporary tag, comm_to_node is now inconsistent

}


void BarberModularity::insert(int node, int comm, double degree_in_comm){

	//D( cout << quality_type << " insert:: node " << node << " into " << comm << " knc= " << degree_in_comm << " kn = " << G.degrees[node] <<
	//" in[" << comm << "] " << in[comm] << " -> " << in[comm] + degree_in_comm << " tot[" << comm + ( (node<G.B)?0:num_communities ) << "] " 
	//<< tot[comm+( (node<G.B)?0:num_communities )] << " -> " << tot[comm+( (node<G.B)?0:num_communities )] + G.degrees[node] << endl );


	in[comm]  += degree_in_comm;
	if( node < G.B ){
		tot[comm] += G.degrees[node];
	} else {
		tot[num_communities + comm] += G.degrees[node];
	}
	node_to_comm[node] = comm; //temporary tag, comm_to_node is now inconsistent


}

double BarberModularity::gain(int node, int comm, double degree_in_comm) {

	//D( cout << quality_type << " gain:: node " << node << " into " << comm << " knc= " << degree_in_comm << " kn = " << G.degrees[node] <<
	//" degree[" << node << "] = " << G.degrees[node]);
	if(node<G.B){
		//D(cout << " tot[" << num_communities+comm << "] = " << tot[num_communities + comm] << " dQ = " << (degree_in_comm - 2*G.degrees[node]*tot[num_communities + comm]/G.num_links) << endl );
		return (degree_in_comm - 2*G.degrees[node]*tot[num_communities + comm]/G.num_links); //factor 2 is for double counting total number of links
	} 
	//D(cout << " tot[" << comm << "] = " << tot[comm] << " dQ = " << (degree_in_comm - 2*tot[comm]*G.degrees[node]/G.num_links) << endl );	
	return (degree_in_comm - 2*tot[comm]*G.degrees[node]/G.num_links); //factor 2 is for double counting total number of links

}

void BarberModularity::join(int comm1, int comm2, double v12){

	D(cout << quality_type << "join " << comm1 << " into " << comm2 << " v12 = " << v12 << endl);

	in[comm2]  += in[comm1] + v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm2+num_communities] += tot[comm1+num_communities];
	tot[comm1] = 0;
	tot[comm1+num_communities] = 0;
	
	BarberModularity::qjoin(comm1, comm2);
  
}





double BarberModularity::join_gain(int comm1, int comm2, double v12) {

	//D( cout << quality_type << " agg_gain:: comm1 " << comm1 << " into " << comm2 << " v12 = " << v12  
	//<< " tot[" << comm1 << "] = " << tot[comm1] << 
	//" tot[" << comm2 << "] = " << tot[comm2] <<
	//" tot[" << comm1+num_communities << "] = " << tot[comm1+num_communities] <<
	//" tot[" << comm2+num_communities << "] = " << tot[comm2] <<
	//" dQ = " << (v12 - (tot[comm1]*tot[comm2+num_communities] + tot[comm1+num_communities]*tot[comm2])*(2/G.num_links)) << endl);

	 //factor 2 is for double counting total number of links
	return (v12 - (tot[comm1]*tot[comm2+num_communities] + tot[comm1+num_communities]*tot[comm2])*(2/G.num_links)); //*(2/G.num_links)

}

#endif
