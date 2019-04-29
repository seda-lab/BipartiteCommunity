// File: graph_binary.h
// -- graph handling header file
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


#ifndef BINARY_GRAPH_H
#define BINARY_GRAPH_H

#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>


using namespace std;


class BinaryGraph {
 public:
  int nb_nodes;
  unsigned long long nb_links;

  long double total_weight;
  int sum_nodes_w;

  vector<unsigned long long> degrees;
  vector<int> links;
  vector<long double> weights; //graphs are now always weighted

  vector<int> nodes_w;

  int bipartite;
  
  BinaryGraph();
  
  // binary file format is
  // 4 bytes for the number of nodes in the graph
  // 8*(nb_nodes) bytes for the cumulative degree for each node:
  //    deg(0)=degrees[0]
  //    deg(k)=degrees[k]-degrees[k-1]
  // 4*(sum_degrees) bytes for the links
  // IF WEIGHTED, 10*(sum_degrees) bytes for the weights in a separate file
  BinaryGraph(char *filename, int bipartite);
  
  vector<vector<pair<int, long double> > > generate_links();
  void generate_from_links(vector<vector<pair<int, long double> > > m, int bipartite_);
  void display_binary(char *outfile);
  void display();

  // assign a weight to a node (needed after the first level)
  void assign_weight(int node, int weight);

  // return the number of neighbors (degree) of the node
  inline int nb_neighbors(int node);

  // return the number of self loops of the node
  inline long double nb_selfloops(int node);

  // return the weighted degree of the node
  inline long double weighted_degree(int node);

  // return pointers to the first neighbor and first weight of the node
  inline pair<vector<int>::iterator, vector<long double>::iterator > neighbors(int node);

};


inline int
BinaryGraph::nb_neighbors(int node) {
  assert(node>=0 && node<nb_nodes);

  if( (node==0) || (node==bipartite) ){
    return degrees[node];
  } else{
    return (int)(degrees[node]-degrees[node-1]);
  }
  
}

inline long double
BinaryGraph::nb_selfloops(int node) {
  assert(node>=0 && node<nb_nodes);

  pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
  for (int i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
		return (long double)*(p.second+i);
    }
  }
  return 0.0L;
}

inline long double
BinaryGraph::weighted_degree(int node) {
  assert(node>=0 && node<nb_nodes);
  

    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    long double res = 0.0L;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      res += (long double)*(p.second+i);
    }
    return res;
  
}

inline pair<vector<int>::iterator, vector<long double>::iterator >
BinaryGraph::neighbors(int node) {
  assert(node>=0 && node<nb_nodes);
  //general idea
  //degrees[node - 1] = number of edges in list, lets us jump to the right place in 'links'
  //when graph is bipartite this is a little more complicated...
  if(bipartite > 0){
	  if(node == 0){
         return make_pair(links.begin(), weights.begin());
      } else if(node == bipartite){
		 return make_pair(links.begin()+degrees[bipartite-1], weights.begin()+degrees[bipartite-1]);
	  } else {
		return make_pair(links.begin() + (node>bipartite?degrees[bipartite-1]:0)+degrees[node-1], weights.begin()+ (node>bipartite?degrees[bipartite-1]:0)+degrees[node-1]);
      } 
  } else {
	    if (node==0){
			return make_pair(links.begin(), weights.begin());
		} else {
			return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
		}
  }
 
}


#endif // BINARY_GRAPH_H
