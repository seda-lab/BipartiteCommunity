#ifndef BLOCKMODULARITY_H
#define BLOCKMODULARITY_H

#include "quality.h"
#include <unordered_set>
#include <cmath>

using namespace std; 

template<typename T>
class BlockModularity : public Quality<T> {
	public:
	
		//This is necessary for a templated class
		//Doing this way avoids having to sue this-> in the method definitions
		using Quality<T>::num_communities;
		using Quality<T>::val;
		using Quality<T>::gamma;
		using Quality<T>::node_to_comm;
		using Quality<T>::comm_to_node;
		using Quality<T>::tot;
		using Quality<T>::in; //the diagonals
		using Quality<T>::dc;
		using Quality<T>::G; 
		using Quality<T>::quality_type;

		//vector< vector<double> > sigma_in;
		unordered_map<int, double> tot2;
		double sq_sum; //sum of squares of rows of Q_ab
		WeightedNetwork GI;
		
		BlockModularity(T& _G);
		BlockModularity(T& _G, unordered_map<int, int> &labels);
	
		void init();
		double eval();
		
		void remove(int node, int comm, double degree_in_comm);
		void insert(int node, int comm, double degree_in_comm);
		void join(int comm1, int comm2, double v12);

		double gain(int node, int comm, double degree_in_comm);
		double join_gain(int comm1, int comm2, double v12);
		
		void print_basic(bool use_original_ids=true);
};


template <typename T>
BlockModularity<T>::BlockModularity(T& _G) : Quality<T>(_G){
	quality_type = "block_modularity";	
	BlockModularity<T>::init();
}

template <typename T>
BlockModularity<T>::BlockModularity(T& _G, unordered_map<int, int> &labels) : Quality<T>(_G, labels){
	quality_type = "block_modularity";	
	BlockModularity<T>::init();
}

///set up in and tot vectors
template <typename T>
void BlockModularity<T>::init(){

	//num_communities += 1;
	/*in.resize( num_communities ); 
	tot.resize( num_communities ); 
	sigma_in.resize( num_communities ); 
	for(int i=0; i<num_communities; ++i){ 
		in[i] = 0; tot[i] = 0;
		sigma_in[i].resize( num_communities ); 
		for(int j=0; j<num_communities; ++j){ sigma_in[i][j] = 0; }		
	}

	for(unsigned i=0;i<G.get_nbrs.size();++i){
		for(unsigned j=0; j<G.get_nbrs[i].size(); ++j){
			int nj = G.get_nbrs[i][j].first;
			sigma_in[ node_to_comm[i] ][ node_to_comm[nj] ] += G.get_nbrs[i][j].second; 
		}
	}

	//loop over nodes
	for(unsigned i=0; i<G.degrees.size(); ++i){
		tot[ node_to_comm[i] ] += G.degrees[i];
	}

	G = G.induced_graph(node_to_comm);
	tot2.resize( G.degrees.size() );
	for(unsigned i=0; i<G.degrees.size(); ++i){
		tot2[ i ] += G.degrees[i];
	}
	dsum = 0;
	for(unsigned i=0; i<G.degrees.size(); ++i){
		in[i] += tot2[ i ]*tot2[ i ];
		dsum += tot2[ i ]*tot2[ i ];
	}*/

	GI = G.induced_graph(node_to_comm);	
	sq_sum = 0;
	for(auto &i : GI.degrees){ tot[i.first] += i.second; }
	for(auto &i : tot){ tot2[i.first] += i.second*i.second; sq_sum += i.second*i.second; }

}


template <typename T>
void BlockModularity<T>::print_basic(bool use_original_ids){

	cout << quality_type << " is the quality function" << endl;
	cout << num_communities << " communities" << endl;
	for(auto &i : node_to_comm){ cout << i.first << " in " << i.second << endl; }
	for(auto &i : comm_to_node){cout << "community " << i.first << " = (";for(auto &j : i.second){cout << j << " ";} cout << ")" << endl;}
	
	for(auto &i : tot){ cout << "tot[" << i.first << "] = " << i.second << " tot2[" << i.first << "] = " << tot2[i.first] << endl; }	
	cout << "sq_sum = " << sq_sum << endl;
	cout << "Current Quality = " << val << endl;
	
}


template <typename T>
double BlockModularity<T>::eval(){
	
	double q = 0;
	for(auto &i : GI.get_nbrs){

		double rowsum = sq_sum;
		for(auto &j : i.second){
			rowsum -= tot2[j.first];
			double add = j.second - (gamma*tot[i.first]*tot[j.first]/GI.num_links);
			q += add*add;		
			//cout << "EVAL::" << "Q_{" << i.first << ";" << j.first << "} = " << j.second << " - " <<  tot[i.first] << "*" << tot[j.first] << "/" << GI.num_links << endl; 				
		}
		//delete after debug
		double tmp = 0;
		for(auto &j : GI.get_nbrs){
			if( i.second.find(j.first) == i.second.end() ){
				//cout << "EVAL::" << "Q_{" << i.first << ";" << j.first << "} = 0 - " <<  tot[i.first] << "*" << tot[j.first] << "/" << GI.num_links << endl; 				
				tmp += pow(gamma*tot[i.first]*tot[j.first]/GI.num_links,2);
				q += pow(gamma*tot[i.first]*tot[j.first]/GI.num_links,2);
			}
		}
		//cout << "EVAL::" << "check " << tmp << " == " << (gamma*gamma)*(tot[i.first]*tot[i.first])*(rowsum)/(GI.num_links*GI.num_links) << endl;
		//q += (gamma*gamma)*(tot[i.first]*tot[i.first])*(rowsum)/(GI.num_links*GI.num_links);
	}
	val = q/(GI.num_links*GI.num_links);
	return val; 
	
}

template <typename T>
void BlockModularity<T>::remove(int node, int comm, double degree_in_comm){

	//remove node from comm and put it on its own at the end
	node_to_comm[node] = -1;
	
	comm_to_node[comm].erase(node);
	comm_to_node[-1] = unordered_set<int>{ node };

	sq_sum -= tot2[comm];
	tot2[comm] -= G.degrees[node]*(2*tot[comm] - G.degrees[node]);
	tot[comm] -= G.degrees[node];
	tot[-1] = G.degrees[node];
	tot2[-1] = G.degrees[node]*G.degrees[node];
	sq_sum += tot2[-1] + tot2[comm];
	
	if( comm_to_node[comm].empty() ){ 
		comm_to_node.erase(comm); 
		tot.erase(comm);
		tot2.erase(comm);
	}

	num_communities = comm_to_node.size();
	
	GI = G.induced_graph( node_to_comm ); //TODO very slow

	
}

template <typename T>
void BlockModularity<T>::insert(int node, int comm, double degree_in_comm){

	//insert isolated node into comm, get rid of extra comm
	node_to_comm[node] = comm;
	
	comm_to_node[comm].insert(node);
	comm_to_node.erase(-1);
	
	sq_sum -= tot2[comm] + tot2[-1];
	tot2[comm] += G.degrees[node]*(2*tot[comm] + G.degrees[node]);
	tot[comm] += G.degrees[node];

	tot.erase(-1);
	tot2.erase(-1);
	sq_sum += tot2[comm];

	num_communities = comm_to_node.size();
	
	GI = G.induced_graph( node_to_comm ); //TODO very slow	
  
}

template <typename T>
double BlockModularity<T>::gain(int node, int comm, double degree_in_comm) {

	return 0;

	/*double q = 0;
	int n = node;
	int c = comm;

	cout << "Gain from inserting " << n << " into " << c << endl;
	for(auto &i : GI.get_nbrs){
		int a = i.first;
		if(a != c){
			cout << "Gain:: " << " checking comm " << a << endl;
			double Qac = 0;
			if( i.second.find(c) == i.second.end()){ //c is not a neighbour of a
				cout << "Q_{" << a << ";" << c << "} = 0 - " <<  tot[a] << "*" << tot[c] << "/" << GI.num_links << endl; 
				Qac = - (gamma*tot[a]*tot[c]/GI.num_links);
			} else {	//c is a neighbour of a
				cout << "Q_{" << a << ";" << c << "} = " << i.second[c] << " - " <<  tot[a] << "*" << tot[c] << "/" << GI.num_links << endl; 				
				Qac = i.second[c] - (gamma*tot[a]*tot[c]/GI.num_links);
			}
			double Qan = 0;
			if( i.second.find(n) == i.second.end()){ //n is not a neighbour of a
				cout << "Q_{" << a << ";" << n << "} = 0 - " <<  tot[a] << "*" << G.degrees[n] << "/" << GI.num_links << endl; 				
				Qan = - (gamma*tot[a]*G.degrees[n]/GI.num_links);
			} else {	//n is a neighbour of a
				cout << "Q_{" << a << ";" << n << "} = " << i.second[n] << " - " <<  tot[a] << "*" << G.degrees[n] << "/" << GI.num_links << endl; 								
				Qan = i.second[n] - (gamma*tot[a]*G.degrees[n]/GI.num_links);
			}
			cout << "Q_{" << a << ";" << c << "} Q_{" << a << ";" << n << "} = " << Qan*Qac << endl;
			q += Qan * Qac;
		}
	}
	q *= 4;
	cout << "EXTRA" << endl;
	//Then add the cn cross terms
	double Qcc = ( GI.get_nbrs[c].find(c) == GI.get_nbrs[c].end() ? 0 : GI.get_nbrs[c][c] ) - (gamma*tot[c]*tot[c]/GI.num_links);
	double Qnn = ( G.get_nbrs[n].find(n) == G.get_nbrs[n].end() ? 0 : G.get_nbrs[n][n] ) - (gamma*G.degrees[n]*G.degrees[n]/GI.num_links);
	double Qcn = 0; for(auto &i : comm_to_node[c]){ if( G.get_nbrs[n].find(i) != G.get_nbrs[n].end() ){Qcn += G.get_nbrs[n][i];} }
	cout << "Q_{" << c << ";" << n << "} = " << Qcn << " - " <<  G.degrees[n] << "*" << tot[c] << "/" << GI.num_links << endl; 		

	Qcn -= (gamma*G.degrees[n]*tot[c]/GI.num_links);
	
	cout << "Q_{" << n << ";" << n << "} = " << ( G.get_nbrs[n].find(n) == G.get_nbrs[n].end() ? 0 : G.get_nbrs[n][n] ) << " - " <<  G.degrees[n] << "*" << G.degrees[n] << "/" << GI.num_links << endl; 		
	cout << "Q_{" << c << ";" << c << "} = " << ( GI.get_nbrs[c].find(c) == GI.get_nbrs[c].end() ? 0 : GI.get_nbrs[c][c] ) << " - " <<  tot[c] << "*" << tot[c] << "/" << GI.num_links << endl; 		
	
	q += (Qnn + Qcc + 2*Qcn)*(Qnn + Qcc + 2*Qcn) - Qnn*Qnn - Qcc*Qcc - 2*Qcn*Qcn;
	
	
	return q; */
}


template <typename T>
void BlockModularity<T>::join(int comm1, int comm2, double v12){

	//D(cout << quality_type << "join " << comm1 << " into " << comm2 << " v12 = " << v12 << endl);

	in[comm2]  += in[comm1] + 2*v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;
	
	BlockModularity<T>::qjoin(comm1, comm2);
  
}


template <typename T>
double BlockModularity<T>::join_gain(int comm1, int comm2, double v12) {

	//D( cout << quality_type << " agg_gain:: comm1 " << comm1 << " into " << comm2 << " v12 = " << v12  
	//<< " tot[" << comm1 << "] = " << tot[comm1] << " tot[" << comm2 << "] = " << tot[comm2] <<
	//" dQ = " << (2*v12 - 2*tot[comm1]*tot[comm2]/G.num_links) );

	return (v12 - gamma*tot[comm1]*tot[comm2]/G.num_links); //*(2/G.num_links)*(2)

}
#endif
