// File: quality.h
// -- quality functions source file
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


#include "quality.h"

using namespace std;


Quality::~Quality() {
  n2c.clear();
}

int
Quality::renumber_comms(){ //renumber communities, return number of communities

  int last=0;
  vector<int> renumber(size, -1);

  // Renumber communities
  for (int node=0 ; node < size ; node++){renumber[n2c[node]]++;}
  for (int i=0 ; i < size ; i++) { if (renumber[i]!=-1){ renumber[i] = last++; } }
  for (int i=0 ; i < size ; i++) { int old_comm = n2c[i]; n2c[i] = renumber[old_comm]; }
  
  return last;
}

int
Quality::fill_c2n(){ //community to node map

  int ncomms = renumber_comms();
  c2n.resize(ncomms);
  for(int c=0; c<ncomms; ++c){ c2n[c].resize(0); }

  for(int node=0; node<size; ++node){ 
	  c2n[ n2c[node] ].push_back(node); 
  }

	
  return ncomms;
}
