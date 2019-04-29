// File: graph.cpp
// -- simple graph handling source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#include "graph.h"
#include <map>
#include <set>
#include <utility> 

using namespace std;


Graph::Graph(char *filename, bool type, int bipartite_) {
	
  ifstream finput;
  finput.open(filename,fstream::in);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  weighted = type;
  bipartite = bipartite_;

  unsigned long long nb_links = 0ULL;

  while (!finput.eof()) {
    unsigned int src, dest;
    long double weight = 1.0L;

    if (weighted) {
      finput >> src >> dest >> weight;
    } else {
      finput >> src >> dest;
    }
	
    if (finput) {

      if (links.size()<=max(src,dest)+1) {
        links.resize(max(src,dest)+1);
      }
      
      links[src].push_back(make_pair(dest,weight));
		
      if(src!=dest){
        links[dest].push_back(make_pair(src,weight));
      } 
      nb_links += 1ULL;
    }

  }
  finput.close();


}


void
Graph::clean() {
  for (unsigned int i=0 ; i<links.size() ; i++) {
    map<int, long double> m;	
    map<int, long double>::iterator it;

    for (unsigned int j=0 ; j<links[i].size() ; j++) { //for every src
      it = m.find(links[i][j].first);	//search m for dest
      if (it==m.end()){					//if not found insert
		  m.insert(make_pair(links[i][j].first, links[i][j].second));
	  } else {							//if found add weight to link
		weighted=true;
      	it->second += links[i][j].second;
	  }
    }
    
    vector<pair<int, long double> > v; //new links vector
    for (it = m.begin() ; it!=m.end() ; it++){ v.push_back(*it); }
    links[i].clear();
    links[i] = v;
  }
}

bool
Graph::renumber(char *filename) {

  set<int> uniq; //number of unique nodes
  int min_idx = -1;
  for(unsigned int i=0; i<links.size(); ++i){ 
	  if( (min_idx<0) && (links[i].size()>0) ){min_idx = i;} //first non empty node list
	  for (unsigned int j=0 ; j<links[i].size() ; j++) {
		uniq.insert(links[i][j].first); //insert nodes an end of links
      } 
  }
  
  if( uniq.size() < links.size() && min_idx >= 0){ //fewer unique nodes than ids
	  
	  vector<int> linked(links.size(),-1); //this will = 1 if the id has at least 1 edge
	  renum.resize(links.size()); //this will be the map renum[ new_node_id ] = old_node_id;

	  for(unsigned int i=0; i<renum.size(); ++i){ renum[i] = -1; }
	  int nb = 0;

	  ofstream foutput;
	  if( filename ){ foutput.open(filename, fstream::out); }
	  
	  for (unsigned int i=0 ; i<links.size() ; i++) { if (links[i].size() > 0){linked[i] = 1;} }
	  
	  //assign the new ids
	  for (unsigned int i=0 ; i<links.size() ; i++) {
		if (linked[i]==1) { 
		  renum[i] = nb++;		//renum[ old_node_id ] = new_node_id or -1
		  if(filename){ foutput << i << " " << renum[i] << endl; }
		}
	  }

	  //update the links vector
	  for (unsigned int i=0 ; i<links.size() ; i++) {
		if (linked[i]==1) { //has an edge
		  for (unsigned int j=0 ; j<links[i].size() ; j++) {
			links[i][j].first = renum[links[i][j].first];	//change dest ids
		  }
		  links[renum[i]] = links[i]; //change src ids
		}
	  }
	  links.resize(nb); //get rid of empty ids

	  //update the bipartite index
	  if( bipartite > 0 ){ 
		  int new_bip = 0;
		  for(unsigned int i=0; i<renum.size(); ++i){
			  if(renum[i] > -1){
				   if((int)i >= bipartite){ break; }
				   new_bip += 1;
			  }
		  }
		  bipartite = new_bip;
	  }
	  	  
	  vector<int> tmp(nb); //tmp[ new_node_id ] = old_node_id
	  for(unsigned int i=0; i<renum.size(); ++i){ 
		  if( renum[i] != -1){tmp[ renum[i] ] = i;}		  
	  };
	  renum.resize(nb); renum = tmp; //renum[ new_node_id ] = old_node_id

	  return true;
  }
  return false;
  
}

void
Graph::project(int side){

  int from, to;
  if(side == 1){from = 0; to = bipartite; } //project out left nodes
  else{from=bipartite; to=(int)links.size(); }	//project out right nodes

  int max_idx = 0;
  map< pair<int, int>, int > A;  //count the number of 1/0 nodes linking 0/1 nodes
  //computing A_ij = sum_k links_ik links_jk
  for (int i=from ; i<to ; i++) { 
	  
	  for (unsigned int j1=0 ; j1<links[i].size() ; j1++) {		
		  int id1 = links[i][j1].first;
	  for (unsigned int j2=0 ; j2<links[i].size() ; j2++) {	
		  int id2 = links[i][j2].first;	
		  
			pair<int, int> idx = make_pair(id1, id2);
			if( A.find(idx) != A.end() ){
			  A[idx] += links[i][j1].second*links[i][j2].second;
			} else {
			  A[idx] = links[i][j1].second*links[i][j2].second;
			}
			int m = std::max(id1, id2);
		    if( m>max_idx ){ max_idx=m; }		
		      
	}}
  }
  //transform A from map to vector
  vector<vector<pair<int, long double> > > links_(max_idx+1, vector<pair<int, long double> >(0) );
  for(std::map<pair<int, int>, int >::iterator it1=A.begin(); it1!=A.end(); ++it1){
	  if((it1->first).first < (it1->first).second){
		links_[(it1->first).first].push_back( make_pair((it1->first).second, (it1->second)) );
		links_[(it1->first).second].push_back( make_pair((it1->first).first, (it1->second)) );
	}
  }
  
  //After projection it is no longer a bipartite graph!
  set_from(links_, true, 0);
  /*
  //if after projecting we have nodes without edges
  if( renum.size()>0 && naivesize > (int)renum.size() ){
	cerr << "naive size " << naivesize << " actual size " << renum.size() << endl;
	
	vector<int> tmp(naivesize, -1);
	int sub = ((side == 0)?0:bipartite_original);  
	for(unsigned int i=0; i<renum.size();++i){ tmp[ renum[i]-sub ] = 1; } //tmp[node] = 1 if node has an edge
	missing.resize(0);
	for(unsigned int i=0; i<tmp.size(); ++i){ if(tmp[i] < 0){missing.push_back(sub+i); }  } //missing = list of nodes without edges
	//for(unsigned int i=0; i<missing.size(); ++i){ cerr << "missing " << missing[i] << endl; }
  }*/
  projected = side;
}





void
Graph::display(bool ordered) {
  for (unsigned int i=0 ; i<links.size() ; i++) {
    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      int dest = links[i][j].first;
      if( ordered && ((int)i>dest) ){ continue; }
      long double weight = links[i][j].second;
      if (weighted){
		cout << i << " " << dest << " " << weight << endl;
      } else {
		cout << i << " " << dest << endl;
	  }
    }
  }
}

void
Graph::display_txt(bool ordered, char *filename) {
  ofstream foutput;
  foutput.open(filename, fstream::out);

  int sub = 0;
  if( projected >= 0){ sub = ((projected == 0)?0:bipartite); }

  for (unsigned int i=0 ; i<links.size() ; i++) {
	for (unsigned int j=0 ; j<links[i].size() ; j++) {
	  int dest = links[i][j].first;
	  if( ordered && ((int)i>dest) ){ continue; }
	  long double weight = links[i][j].second;
	  if (weighted){
		if( renum.size() > 0 ){
		  foutput << renum[i]-sub << " " << renum[dest]-sub << " " << weight << endl;
		} else {
		  foutput << i << " " << dest << " " << weight << endl;
		}
	  } else {
		if( renum.size() > 0 ){
		  foutput << renum[i]-sub << " " << renum[dest]-sub << endl;
		} else {
		  foutput << i << " " << dest << endl;
		}
	  }
	}
  }
  
  foutput.close();

}

void
Graph::display_binary(char *filename) {
  ofstream foutput;
  foutput.open(filename, fstream::out | fstream::binary);

  int s = links.size();

  // outputs number of nodes
  foutput.write((char *)(&s),sizeof(int));
  
  // outputs cumulative degree sequence
  unsigned long long tot = 0ULL;
  for (int i=0 ; i<s ; i++) {
    tot += (unsigned long long)links[i].size();
    foutput.write((char *)(&tot),sizeof(unsigned long long));
  }

  // outputs links
  for (int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links[i].size() ; j++) {
      int dest = links[i][j].first;
      foutput.write((char *)(&dest),sizeof(int));
    }
  }

    for (int i=0 ; i<s ; i++) {
      for (unsigned int j=0 ; j<links[i].size() ; j++) {
		long double weight = links[i][j].second;
		foutput.write((char *)(&weight),sizeof(long double));
      }
    }
  
  foutput.close();

}
