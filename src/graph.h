// File: graph.h
// -- simple graph handling header file
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


#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;


class Graph {
 public:
  vector<vector<pair<int, long double> > > links;
  vector<int> renum;
  vector<int> missing;
  int bipartite;
  bool weighted;
  int projected;
  
  Graph (){bipartite = 0; projected=-1; weighted=false; renum.resize(0); links.resize(0); }
  Graph (char *filename, bool type, int bipartite);

  void clean(); //amalgamates repeated links
  bool renumber(char *filename); //number nodes 0,1,...,nb_nodes-1. return if renumbered
  void project(int num_nodes); //project onto node set.  return if renumbered
  void display(bool ordered);
  void display_binary(char *filename);
  void display_txt(bool ordered, char *filename);

  void set_from(vector<vector<pair<int, long double> > > links_, bool weighted_, int bipartite_) {
	weighted = weighted_;
	bipartite = bipartite_;
	links = links_;
  }
};

#endif // GRAPH_H
