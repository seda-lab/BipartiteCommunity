// File: graph_binary.cpp
// -- graph handling source
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


#include <fstream>
#include "graph_binary.h"


BinaryGraph::BinaryGraph() {
  nb_nodes = 0;
  nb_links = 0ULL;

  total_weight = 0.0L;
  sum_nodes_w = 0;
  
  bipartite=0;

}

BinaryGraph::BinaryGraph(char *filename, int bipartite_) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  // Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, sizeof(int));
  if (finput.rdstate() != ios::goodbit) {
    cerr << "The file " << filename << " is not a valid graph" << endl;
    exit(EXIT_FAILURE);
  }
  bipartite = bipartite_;

  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);

  finput.read((char *)&degrees[0], nb_nodes*sizeof(unsigned long long));
  if(bipartite>0){
	for(int i=bipartite; i<nb_nodes; i++){
	  degrees[i] -= degrees[bipartite-1];
	}
  }

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links = degrees[nb_nodes-1];
  if(bipartite>0){nb_links *= 2;}
  
  links.resize(nb_links);
  finput.read((char *)(&links[0]), nb_links*sizeof(int));

  // IF WEIGHTED, read weights: 10 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight = 0.0L;
  
  weights.resize(nb_links);
  finput.read((char *)(&weights[0]), nb_links*sizeof(long double));

  // Compute total weight
  for (int i=bipartite ; i<nb_nodes ; i++){ total_weight += (long double)weighted_degree(i); }
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}

void
BinaryGraph::generate_from_links(vector<vector<pair<int, long double> > > m, int bipartite_){
  
  nb_nodes = m.size();
  bipartite = bipartite_;
  
  // outputs cumulative degree sequence
  degrees.resize(nb_nodes);
  unsigned long long tot = 0ULL;
  for (int i=0 ; i<nb_nodes; i++) {
	if(i==bipartite){tot=0;}
    tot += (unsigned long long)m[i].size();
    degrees[i] = tot;
  }

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links = degrees[nb_nodes-1];
  if(bipartite>0){nb_links *= 2;}
 
  links.resize(nb_links);  
  int ct = 0;
  for (int i=0 ; i<nb_nodes ; i++) {
    for (unsigned int j=0 ; j<m[i].size() ; j++) {
      links[ct] = m[i][j].first;
      ct += 1;
    }
  }
  
  weights.resize(nb_links);
  ct = 0;
  for (int i=0 ; i<nb_nodes ; i++) {
    for (unsigned int j=0 ; j<m[i].size() ; j++) {
      weights[ct++] = m[i][j].second;
    }
  }
  
  // Compute total weight
  for (int i=bipartite ; i<nb_nodes ; i++){ 
	  total_weight += (long double)weighted_degree(i); 
  }
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;

}

vector<vector<pair<int, long double> > > 
BinaryGraph::generate_links() {

  vector<vector<pair<int, long double> > > links_(nb_nodes);

  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node <= *(p.first+i)) {
		links_[node].push_back( make_pair( *(p.first+i), *(p.second+i) ) );
		if(node != *(p.first+i)){
			links_[*(p.first+i)].push_back( make_pair( node, *(p.second+i) ) );
		}
      }
    }
  }
  return links_;
}

void
BinaryGraph::assign_weight(int node, int weight) {
  sum_nodes_w -= nodes_w[node];

  nodes_w[node] = weight;

  sum_nodes_w += weight;
}

void
BinaryGraph::display() {
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    cout << node << ":" ;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (true) {
	if (weights.size()!=0)
	  cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
	else
	  cout << " " << *(p.first+i);
      }
    }
    cout << endl;
  }
}

void
BinaryGraph::display_binary(char *outfile) {
  ofstream foutput;
  foutput.open(outfile ,fstream::out | fstream::binary);

  foutput.write((char *)(&nb_nodes),sizeof(int));
  foutput.write((char *)(&degrees[0]),sizeof(unsigned long long)*nb_nodes);
  foutput.write((char *)(&links[0]),sizeof(int)*nb_links);
  foutput.write((char *)(&weights[0]),sizeof(long double)*nb_links);
}

